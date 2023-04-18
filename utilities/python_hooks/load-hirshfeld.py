"""
This aims Python hook replaces the Hirshfeld volumes with those loaded from a
previous aims output (in hirshfeld.out). It needs to be run in parallel.
"""
import numpy as np

volumes = None


# this is called when parsing control.in
def parse():
    global volumes  # otherwise volumes would be local to function
    with open('hirshfeld.out') as f:
        while 'Performing Hirshfeld analysis' not in next(f):
            pass  # find the Hirshfeld section
        volumes = np.array(list(get_volumes(f)))


def get_volumes(f):
    for line in f:
        if not line.strip():
            return
        if 'Free atom volume' in line:
            free = float(line.split()[-1])
        elif 'Hirshfeld volume' in line:
            hirsh = float(line.split()[-1])
            yield hirsh/free


# this will be called right after Hirshfeld analysis.
def run(ctx):
    assert volumes.size == ctx.coords.shape[1]  # check we have correct shape
    ctx.hirshfeld_volume[:] = volumes  # replace values in-place
