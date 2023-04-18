#!/usr/bin/python

import re
import difflib
import os
import os.path
import shutil
import subprocess
import optparse
import tempfile
import sys
import array

ABSTOL = 1e-10
RELTOL = 1e-7

FLOAT_RE_STR = r'''
  [-+]?      # optional sign
  (?:        # either:
    \d+\.\d*(?:[eE][-+]?\d+)? | # 1.e1 | 1.
    \d+[eE][-+]?\d+           | # 1e1
    \.\d+(?:[eE][-+]?\d+)?      # .1e1 | .1
  )
'''
FLOAT_RE = re.compile(FLOAT_RE_STR, re.X)
COMPLETE_FLOAT_RE = re.compile("^" + FLOAT_RE_STR + "$", re.X)

USAGE = """
  %prog [options] file1 file2
  %prog [options] dir1 dir2

Prepare two files or two directories for comparison.  Thereby
ignore all differences in given regular expressions as well
as numerical differences smaller than atol + rtol*(|x|+|y|)/2
where x and y are to be compared.

The config file syntax is simple and line based.  Lines starting
with '#' are comments, the text after an initial 'ignore:' should be
a regular expression.

This script is experimental.  Please mail any comments to
wieferink@fhi-berlin.mpg.de.
"""

class Matches(object):
    """Keep book of ranges within a string where some pattern matches.

    Example:
    >>> text = 'This foo text is bar.'
    >>> obj = Matches([r'foo', r'bar'], text)
    >>> obj.contains((4, 4)) # After "This"
    False
    >>> obj.contains((7, 7)) # Within "foo"
    (5, 8)
    >>> text[3:7]
    's fo'
    >>> obj.contains((3, 7))
    False
    >>> obj.contains((17, 19)) # 'ba' within 'bar'
    (17, 20)
    """
    def __init__(self, patterns, text):
        self.intervals = []
        for pat in patterns:
            if isinstance(pat, (str, unicode)):
                pat = re.compile(pat)
            for match in pat.finditer(text):
                self.intervals.append((match.start(), match.end(), pat))
    def contains(self, interval):
        for start, end, pat in self.intervals:
            if interval[0] >= start and interval[1] <= end:
                return (start, end) # which is True
        return False

class NumDiffer(object):
    """Environment object for numerical accuracy aware differ.

    The main entry method is diff_nodes, which can compare either directories
    or files."""
    def __init__(self, atol=ABSTOL, rtol=RELTOL, ignore_patterns=None,
                 shortcut_dirs=False, marker_pattern=None):
        self.atol = atol
        self.rtol = rtol
        if ignore_patterns is None:
            self.ign = []
        else:
            self.ign = ignore_patterns
        self.mark = marker_pattern
        self.shortcut_dirs = shortcut_dirs
        self.numdiffs = dict()
        self.current_file = None
        self.numdiffs[None] = []
        self.n_linediffs = 0

    def tol(self, mag):
        return mag*self.rtol + self.atol

    def _numbers_differ(self, x, y):
        mag = 0.5 * (abs(x) + abs(y))
        err = abs(x - y)
        if err > 0.:
            self.numdiffs[self.current_file].append((mag, err))
        # print x, y, err > self.tol(mag)
        return err > self.tol(mag)

    def tokenize_string(self, line):
        """Tokenize line in ignored, number and other sections.

        >>> pat = re.compile(r'aa.zz')
        >>> differ = NumDiffer(ignore_patterns=[pat])
        >>> assert differ.tokenize_string('foo aaizzf 55.') == [
        ...    ('text', 'foo ', 0, 4, 'foo '),
        ...    ('ignore', pat, 4, 9, 'aaizz'),
        ...    ('text', 'f ', 9, 11, 'f '),
        ...    ('number', None, 11, 14, '55.')]
        """

        # === tokenize ignores
        igns = Matches(self.ign, line)
        # sort to make them greedy: (1, 19) < (1, 15) < (2, 50)
        sorted_intervals = sorted(igns.intervals, key=lambda x:(x[0], -x[1]))
        pos = 0
        res = []
        for start, stop, pattern in sorted_intervals:
            if start > pos:
                text = line[pos:start]
                res.append(('text', text, pos, start, text))
                pos = start
            if start == pos:
                res.append(('ignore', pattern, start, stop, line[pos:stop]))
                pos = stop
        if pos < len(line):
            text = line[pos:]
            res.append(('text', text, pos, len(line), text))
        
        # === tokenize numbers
        nums = Matches([FLOAT_RE], line)
        i = 0
        while i < len(res):
            tag, value, tstart, tstop, text = res[i]
            if tag != 'text':
                i += 1
                continue
            assert value == text
            pos = tstart
            fres = []
            for m in FLOAT_RE.finditer(text):
                numstr = text[m.start():m.end()]
                fstart = m.start() + tstart
                fstop = m.end() + tstart
                if fstart > pos:
                    text = line[pos:fstart]
                    fres.append(('text', text, pos, fstart, text))
                    pos = fstart
                numtext = line[fstart:fstop]
                fres.append(('number', None, pos, fstop, numtext))
                pos = fstop
            if pos < tstop:
                text = line[pos:tstop]
                fres.append(('text', text, pos, tstop, text))
                pos = tstop
            res[i:i+1] = fres
            i += len(fres)
        return res
        


    def diff_strings(self, line1, line2):
        """Compare two lines for significant differences.

        >>> differ = NumDiffer(ignore_patterns=[re.compile(r'aa.zz')])
        >>> line1 = 'Test string; x=1.000000000001, aaizz'
        >>> line2 = 'Test string; x=1.000000000002, aajzz'
        >>> differ.diff_strings(line1, line2)
        False
        >>> differ.diff_strings('a', 'b')
        True
        """
        tokens1 = self.tokenize_string(line1)
        tokens2 = self.tokenize_string(line2)
        diff_tokens1 = [(tag, value)
                        for tag, value, start, stop, text in tokens1]
        diff_tokens2 = [(tag, value)
                        for tag, value, start, stop, text in tokens2]
        if diff_tokens1 != diff_tokens2:
            return True # significant difference
        actually_differ = False
        for tok1, tok2 in zip(tokens1, tokens2):
            tag1, value1, start1, stop1, text1 = tok1
            tag2, value2, start2, stop2, text2 = tok2
            if tag1 == 'number':
                if self._numbers_differ(float(text1), float(text2)):
                    actually_differ = True
        return actually_differ

            
    def _diff_lines(self, line1, line2):
        if self.diff_strings(line1, line2):
            res1 = ["#" + line1]
            res2 = ["#" + line2]
            self.n_linediffs += 1
        else:
            res1 = res2 = ["-" + line1, "+" + line2]
        return res1, res2

    def _match_differing_linelists(self, linelist1, linelist2):
        """Compare lists of lines which are expected to differ."""
        left1, left2 = [], []
        right1, right2 = [], []
        differ = difflib.SequenceMatcher(None)
        while True:
            # End of recursion (empty list)?
            if len(linelist1) == 0 or len(linelist2) == 0:
                self.n_linediffs += len(linelist1) + len(linelist2)
                return (left1 + ["|" + line for line in linelist1] + right1,
                        left2 + ["|" + line for line in linelist2] + right2)
            # compare first lines
            differ.set_seqs(linelist1[0], linelist2[0])
            if differ.ratio() > 0.5:
                line1, line2 = linelist1.pop(0), linelist2.pop(0)
                res1, res2 = self._diff_lines(line1, line2)
                left1.extend(res1)
                left2.extend(res2)
                continue
            # compare last lines
            differ.set_seqs(linelist1[-1], linelist2[-1])
            if differ.ratio() > 0.5:
                line1, line2 = linelist1.pop(), linelist2.pop()
                res1, res2 = self._diff_lines(line1, line2)
                right1[0:0] = res1
                right2[0:0] = res2
                continue
            break
        
        # The difficult case: Find best match.
        ratios = []
        for i2, line2 in enumerate(linelist2):
            differ.set_seq2(line2)
            for i1, line1 in enumerate(linelist1):
                differ.set_seq1(line1)
                ratios.append(differ.ratio())
        maxratio = max(ratios)
        max_n = ratios.index(maxratio)
        max_i1 = max_n % len(linelist1)
        max_i2 = max_n / len(linelist1)
        
        # recurse
        l1, l2 = self._match_differing_linelists(linelist1[:max_i1],
                                                 linelist2[:max_i2])
        c1, c2 = self._diff_lines(linelist1[max_i1], linelist2[max_i2])
        r1, r2 = self._match_differing_linelists(linelist1[max_i1+1:],
                                                 linelist2[max_i2+1:])
        return (left1 + l1 + c1 + r1 + right1,
                left2 + l2 + c2 + r2 + right2)

    def _diff_linelists(self, linelist1, linelist2):
        """Compare general lists of lines."""
        differ = difflib.SequenceMatcher(None, linelist1, linelist2)
        codes = differ.get_opcodes()

        reslist1 = []
        reslist2 = []
        for tag, start1, stop1, start2, stop2 in codes:
            lines1 = linelist1[start1:stop1]
            lines2 = linelist2[start2:stop2]
            if tag == 'equal':
                reslist1.extend([" " + line for line in lines1])
                reslist2.extend([" " + line for line in lines2])
            elif tag == 'delete':
                reslist1.extend(["<" + line for line in lines1])
                self.n_linediffs += len(lines1)
            elif tag == 'insert':
                reslist2.extend([">" + line for line in lines2])
                self.n_linediffs += len(lines2)
            elif tag == 'replace':
                res1, res2 = self._match_differing_linelists(lines1, lines2)
                reslist1.extend(res1)
                reslist2.extend(res2)
            else:
                raise ValueError("Invalid opcode %s" % tag)
        return reslist1, reslist2

    def _diff_marked_linelists(self, linelist1, linelist2):
        """Compare lists of lines after searching for markers."""
        reslist1 = []
        reslist2 = []
        while linelist1 or linelist2:
            acclist1 = []
            acclist2 = []
            while linelist1:
                line1 = linelist1.pop(0)
                acclist1.append(line1)
                if self.mark.match(line1):
                    break
                line1 = ""
            while linelist2:
                line2 = linelist2.pop(0)
                acclist2.append(line2)
                if self.mark.match(line2):
                    break
                line2 = ""
            if line1 != line2:
                break
            r1, r2 = self._diff_linelists(acclist1, acclist2)
            reslist1.extend(r1)
            reslist2.extend(r2)
            acclist1 = []
            acclist2 = []
        acclist1.extend(linelist1)
        acclist2.extend(linelist2)
        if acclist1 or acclist2:
            r1, r2 = self._diff_linelists(acclist1, acclist2)
            reslist1.extend(r1)
            reslist2.extend(r2)
        return reslist1, reslist2


    def diff_nodes(self, path1, path2, topath1, topath2, recursive=False):
        """Diff dirs/files path1, path2 (recursively)."""
        if os.path.isdir(path1) and os.path.isdir(path2):
            os.mkdir(topath1)
            os.mkdir(topath2)
            subnodes1 = set(os.listdir(path1))
            subnodes2 = set(os.listdir(path2))
            common = subnodes1 & subnodes2
            # simply copy files without partner
            for subnode in subnodes1 - common:
                subprocess.call("cp -r %s/%s %s" % (path1, subnode, topath1),
                                shell=True)
            for subnode in subnodes2 - common:
                subprocess.call("cp -r %s/%s %s" % (path2, subnode, topath2),
                                shell=True)
            for subnode in common:
                self.diff_nodes(os.path.join(path1, subnode),
                                os.path.join(path2, subnode),
                                os.path.join(topath1, subnode),
                                os.path.join(topath2, subnode),
                                recursive=True)
        elif self.shortcut_dirs and recursive:
            subprocess.call("cp -r %s %s" % (path1, topath1), shell=True)
            subprocess.call("cp -r %s %s" % (path2, topath2), shell=True)
        else:
            self.current_file = os.path.basename(path1)
            if self.current_file not in self.numdiffs:
                self.numdiffs[self.current_file] = []
            file1 = open(path1)
            file2 = open(path2)
            lines1 = file1.readlines()
            lines2 = file2.readlines()
            file1.close()
            file2.close()
            if self.mark is not None:
                lines1, lines2 = self._diff_marked_linelists(lines1, lines2)
            else:
                lines1, lines2 = self._diff_linelists(lines1, lines2)
            file1 = open(topath1, "w")
            file2 = open(topath2, "w")
            file1.writelines(lines1)
            file2.writelines(lines2)
            file1.close()
            file2.close()


def _test():
    import doctest
    try:
        doctest.testmod(raise_on_error=True)
    except doctest.DocTestFailure, failure:
        print 'DocTestFailure:'
        print '     source:', failure.example.source
        print '   expected:', failure.example.want
        print '     actual:', failure.got
        sys.exit(1)
    except doctest.UnexpectedException, failure:
        print 'UnexpectedException:'
        print '     source:', failure.example.source
        e = failure.exc_info
        raise e[0], e[1], e[2]
        sys.exit(1)

def _main():
    formatter = optparse.IndentedHelpFormatter(max_help_position=32)
    parser = optparse.OptionParser(usage=USAGE, formatter=formatter)
    parser.add_option("-d", "--differ", dest="differ", metavar="PROGRAM",
                      default="diff -ur",
                      help="Programe with which to do the diffs")
    parser.add_option("-a", "--abstol", dest="atol", metavar="atol",
                      action="store", type="float",
                      default=1e-10, help="Absolute tolerance [%g]" % ABSTOL)
    parser.add_option("-r", "--reltol", dest="rtol", metavar="rtol",
                      action="store", type="float",
                      default=1e-10, help="Relative tolerance [%g]" % RELTOL)
    parser.add_option("-i", "--ignore", dest="ignore",
                      action="append", default=[],
                      help="Patterns not to diff")
    parser.add_option("-c", "--config", dest="config",
                      action="append", default=[],
                      help="Config file")
    parser.add_option("--min-diff", action="store", type="int", default=0,
                      help="Minimum number of line diffs to start differ")
    parser.add_option("-p", "--plot", dest="plot", action="store_true",
                      help="Plot differing numbers")
    parser.add_option("-s", "--shortcut-dirs", dest="shortcut_dirs",
                      action="store_true",
                      help="Do not parse files in subdirectories")
    parser.add_option("-t", "--test", dest="test", action="store_true",
                      help="Run unit tests")

    (options, args) = parser.parse_args()

    if options.test:
        _test()
        sys.exit(0)

    try:
        path1, path2 = args
    except ValueError:
        parser.error("Need exactly two arguments")

    ignore_list = options.ignore[:]
    marker_pattern = None
    for fn in options.config:
        f = open(fn)
        for line in f:
            line = line.strip()
            if line.strip().startswith('#'):
                continue
            cmd, sep, value = line.partition(':')
            if cmd == 'ignore':
                ignore_list.append(value.strip())
            elif cmd == 'marker':
                if marker_pattern is not None:
                    parser.error('Only one marker allowed')
                marker_pattern = re.compile(value.strip())
            else:
                parser.error('Unknown command "%s" in %s' % (cmd, fn))
                
    ignore_patterns = [re.compile(ignstr) for ignstr in ignore_list]

    numdiffer = NumDiffer(options.atol, options.rtol, ignore_patterns,
                          options.shortcut_dirs, marker_pattern)

    tempdir1 = tempfile.mkdtemp()
    tempdir2 = tempfile.mkdtemp()
    p = None
    try:
        to1 = os.path.join(tempdir1, os.path.basename(path1.rstrip('/')))
        to2 = os.path.join(tempdir2, os.path.basename(path2.rstrip('/')))
        numdiffer.diff_nodes(path1, path2, to1, to2)
        if numdiffer.n_linediffs >= options.min_diff:
            cmd = "%s %s %s" % (options.differ, to1, to2)
        else:
            cmd = ("echo 'Have n_linediffs < %i; do not launch differ'" %
                   options.min_diff)
        p = subprocess.Popen(cmd, shell=True)

        numdiffs = dict()
        n_diff = 0
        alldiffs = []
        for f in numdiffer.numdiffs:
            numdiffs[f] = numdiffer.numdiffs[f]
            alldiffs.extend(numdiffer.numdiffs[f])
            n_diff += len(numdiffs[f])

        sys.stderr.write("Found %i files with %i differing numbers "
                         "and %i differing lines.\n" %
                         (len(numdiffer.numdiffs), n_diff,
                          numdiffer.n_linediffs))
        if options.plot:
            try:
                import warnings
                warnings.filterwarnings('ignore', category=DeprecationWarning)
                import numpy as np
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                ax.set_xscale('log')
                ax.set_yscale('log')
                if alldiffs:
                    alldiffs = np.asarray(alldiffs)
                    alldiffs[:,0] = np.maximum(alldiffs[:,0], 1e-15)
                    min_v = min(np.min(alldiffs[:,0]), 1e-5)
                    max_v = max(np.max(alldiffs[:,0]), 1e5)
                    tt = np.exp(np.linspace(np.log(min_v), np.log(max_v), 200))
                    ax.plot(tt, options.atol + tt*options.rtol, '-',
                            label="tolerance")
                    ax.plot(tt, tt, '--', color='grey', label='identity')
                    ax.set_xlabel("Magnitude")
                    ax.set_ylabel("Error")
                    #for f in numdiffs:
                    #    if nummdiffs[f]:
                    #        this = np.asarray(numdiffs[f])
                    #        ax.plot(this[:,0], this[:,1],
                    #                ls='', marker='o', label=f)
                    #ax.legend()
                    ax.plot(alldiffs[:,0], alldiffs[:,1], 'ro')
                    ax.set_ylim(ymin=1e-15)
                    def on_q_exit(event):
                        if event.key == "q": plt.close(event.canvas.figure)
                    fig.canvas.mpl_connect('key_press_event', on_q_exit)
                    plt.show()

            except ImportError:
                sys.stderr.write(
                    "With numpy & matplotlib this could be visualized.\n")
        p.wait()
    finally:
        if p is not None:
            p.wait()
        shutil.rmtree(tempdir1)
        shutil.rmtree(tempdir2)


if __name__ == "__main__":
    _main()


