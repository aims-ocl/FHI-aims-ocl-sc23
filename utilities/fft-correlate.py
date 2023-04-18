#!/usr/bin/python

# Reads an output file of FHI-aims, gets the dipole moments and calculates
# various time correlation functions -- MR 2009

import numpy as np
import sys
from optparse import OptionParser


C_LIGHT = 2.99792458e8 # m/s
C_LIGHT_CM_PS = C_LIGHT / 1e-2 * 1e-12

def testmod():
   import doctest
   try:
      doctest.testmod(raise_on_error=True)
   except doctest.DocTestFailure, failure:
      print 'DocTestFailure:'
      print '     source:', failure.example.source.rstrip()
      print '   expected:', failure.example.want.rstrip()
      print '     actual:', failure.got.rstrip()
      sys.exit(1)
   except doctest.UnexpectedException, failure:
      print 'UnexpectedException:'
      print '     source:', failure.example.source.rstrip()
      e = failure.exc_info
      raise e[0], e[1], e[2]
   sys.exit(1)

def error(msg):
   sys.stderr.write(msg + "\n")
   sys.exit(2)

def extend_to_whole(f):
   """Extend array f to one with f_ext[i]==f_ext[2*N-i]==f[i] and f[N]==0.

   The array f is assumed to be the positive part of an even function.
   The resulting array f_ext has the properties:

   >>> from numpy.random import random
   >>> f = random(5)
   >>> f_ext = extend_to_whole(f)
   >>> np.all(f[i] == f_ext[i] == f_ext[-i] for i in range(1, len(f)))
   True
   >>> f_ext[len(f)] == 0.
   True
   """
   N = len(f)
   f_ext = np.zeros(2*N)
   f_ext[0:N] = f
   f_ext[N] = 0.
   f_ext[N+1:2*N] = f[N:0:-1] # all but the zero-th component in reverse order.
   return f_ext

def extend_to_sym(f):
   """Extend array f to one with f_ext[N-1+i]==f_ext[N-1-i]==f[i].

   >>> from numpy.random import random
   >>> f = random(5)
   >>> f_ext = extend_to_sym(f)
   >>> len(f_ext) == 2*len(f) - 1
   True
   >>> f[0] == f_ext[len(f)-1]
   True
   >>> np.all(f[i] == f_ext[len(f)-1+i] == f_ext[len(f)-1-i]
   ...        for i in range(len(f)))
   True
   """
   N = len(f)
   f_ext = np.empty(2*N-1)
   f_ext[:N-1] = f[:0:-1]
   f_ext[N-1] = f[0]
   f_ext[N:] = f[1:]
   return f_ext

def to_high_two(n, qual=3.):
   """Get nearest integer below n which is quite a high power of 2.

   The paramter `qual` adjusts the trade-off between efficiency and accuracy.
   Higher values give lower powers in two, which are nearer to n.

   >>> to_high_two(11)
   8
   >>> to_high_two(13)
   12
   """
   max2 = int(np.log(float(n)) / np.log(2.)) + 1
   best_m = n
   best_q = qual + 1.
   for twopower in xrange(max2):
      fac = 2**twopower
      m = (n // fac) * fac
      approx_qual = 1. - float(m) / n
      power_qual = 1. - np.log(fac) / np.log(n)
      quality = qual*approx_qual + power_qual
      if quality < best_q:
         best_q = quality
         best_m = m
   return best_m


def to_low_primes(n):
   """Get nearest integer below n which is only a power of 2 and 3.

   Unfortunately, scipy.fft is not optimized for any primes other
   than two, so use to_high_two(n) instead for that purpose.

   >>> to_low_primes(11)
   9
   >>> to_low_primes(13)
   12
   """
   max2 = int(np.log(float(n)) / np.log(2.)) + 1
   max3 = int(np.log(float(n)) / np.log(3.)) + 1
   best = 0
   for i2 in range(max2):
      for i3 in range(max3):
         this = 2**i2 * 3**i3
         if this > n:
            exit
         elif this > best:
            best = this
   return best


def autocorr(F, N, use_scipy=None):
   r"""Return autocorrelation (wrt axis=0) of scalar products (wrt axis=1).

   >>> from numpy.random import random
   >>> a1, a2, b1, b2, c1, c2 = random(6)
   >>> F = np.array([[a1,a2], [b1, b2], [c1,c2]])
   >>> f0 = (a1**2 + a2**2 + b1**2 + b2**2 + c1**2 + c2**2) / 3. # three terms
   >>> f1 = (a1*b1 + a2*b2 + b1*c1 + b2*c2)                 / 2. # two terms
   >>> np.allclose(autocorr(F, 2, True), [f0, f1])
   True
   >>> np.allclose(autocorr(F, 1, True), [f0])
   True
   >>> np.allclose(autocorr(F, 2, True), autocorr(F, 2, False))
   True
   """
   f = np.zeros(N)
   F = np.asarray(F)
   lenF, ndf = F.shape
   use_default = use_scipy is None
   if use_default:
      use_scipy = (lenF == to_high_two(lenF))
   if use_scipy:
      try:
         from scipy.signal import fftconvolve
      except ImportError:
         use_scipy = False
   if use_scipy:
      if use_default:
         sys.stderr.write("| Using scipy.signal.fftconvolve.\n")
      ff = np.zeros((N, ndf))
      n_terms = extend_to_sym(np.arange(lenF, 0, -1, dtype=float))
      for j in range(ndf):
         conv_res = fftconvolve(F[:,j], F[::-1,j], 'full')
         conv_norm = conv_res / n_terms # take mean
         assert len(conv_res) == 2*lenF-1 # need [lenF, lenF+N]
         ff[:,j] = conv_norm[lenF-1:lenF+N-1]
      f = np.sum(ff, axis=1)
   else:
      if use_default:
         sys.stderr.write("| Using np.sum for autocorrelation.\n")
      for i in range(N):
         first_ones = F[:lenF-i, :]
         last_ones  = F[i:lenF, :]
         # sp[j] = dot_product(F[j,:], F[i+j,:])
         sp = np.sum(first_ones * last_ones, axis=1)
         f[i] = np.mean(sp) # = sum(sp) / len(sp)
   return f


def Gaussian(n, sigma):
   """Return Gaussian with width sigma centered in array of (odd) length n.

   >>> np.allclose(Gaussian(3, 1e-5), [0., 1., 0.])
   True
   >>> G = Gaussian(3, 1.)
   >>> b, c = G[0], G[1] # border, center
   >>> np.allclose(b / c, np.exp(-1.**2))
   True
   >>> np.allclose(np.convolve([1,3,0], G, 'same'), [c+3*b, b+3*c, 3*b])
   True
   """
   if n % 2 != 1:
      raise ValueError("n (%i) should be odd and positive" % n)
   x0 = float((n-1)/2)
   x = np.arange(0., float(n))
   Gauss = np.exp(-((x-x0)/sigma)**2)
   return Gauss / np.sum(Gauss)

def Hanning(n):
   """Return array of length n with the right half of a Hanning window.

   Near the origin, the Hanning window is nearly unitiy.
   >>> np.allclose(Hanning(10)[0], 1.)
   True
   
   Near the edge, the Hanning window vanishes
   (Hanning(N)[N] would be zero).
   >>> Hanning(10)[-1] < 2. * (2.*np.pi / (2.*10))**2
   True
   >>> np.allclose(extend_to_whole(Hanning(2)), [1., 0.5, 0., 0.5])
   True
   """
   x = np.arange(0., float(n))
   Hann = 0.5 * (1. + np.cos(2.*np.pi* x / (2.*n)))
   return Hann


def read_aims(in_, want_velocities=False, grep_name=None):
   """Read aims output from stream in_ and output F, Dt

   F is of shape (n_timesteps, ndf) and contains either the dipole moments
   or the velocities.  For dipoles ndf=3, for velocities ndf=3*n_atoms.
   Dt is the time step.
   """
   # read
   dipoles = []
   velocities = []
   if grep_name is not None:
      grep_out = open(grep_name, "w")
   for line in in_:
      if "Number of atoms" in line:
         n_atoms = int(line.split()[-1])
         if grep_name is not None: grep_out.write(line)
      elif "Molecular dynamics time step" in line:
         Dt = float(line.split()[-2]) # always in ps -> 1e-12 s
         if grep_name is not None: grep_out.write(line)
      elif "velocity" in line and not "missing" in line:
         fields = line.split()
         i_x = fields.index("velocity") + 1
         velocities.append(map(float, fields[i_x:i_x+3]))
         if grep_name is not None: grep_out.write(line)
      elif "Total dipole moment" in line:
         dipoles.append(map(float, line.split()[-3:]))
         if grep_name is not None: grep_out.write(line)
   if grep_name is not None:
      grep_out.close()

   # Arrange F
   if want_velocities:
      V = np.array(velocities) # shape: (n_atoms*n_timesteps, 3)
      n_timesteps = V.shape[0] / n_atoms
      if n_timesteps * n_atoms != V.shape[0]:
         error("Number of velocities not divisible by number of atoms")
      F = V.reshape((n_timesteps, 3*n_atoms))
   else:
      F = np.array(dipoles)
   return F, Dt


def _main():

   # === Get options

   usage = """\
usage: auto-correlate.py [options] [<aims.out> [<ft> [<autocorr> [<grep>]]]

   Parse an FHI-aims output file (<aims.out>, default: stdin) and calculate
   the autocorrelation function of the dipole moment or the velocities.

   The Fourier transform of the autocorrelation function (the vibrational
   spectrum) is saved to <ft>, the autocorrelation function itself is saved to
   <autocorr>, if the corresponding names are given.  If no output file is
   named (or --plot is given), try to plot.

   If <grep> is given, write all lines from <aims.out> which are actually
   needed for an auto-correlate.py run.
"""
   parser = OptionParser(usage=usage)
   parser.add_option("-v", "--velocities", action="store_true",
                     help="Use velocities instead of dipole moments.")
   parser.add_option("-c", "--cutoff_ratio", type=float, default=0.3,
                     metavar="CUT",
                     help=("Only use a fraction (1.-CUT) " +
                           " of the autocorrelation function [0.3]."))
   parser.add_option("-b", "--broadening", type=float, default=5.,
                     metavar="SIGMA",
                     help="Final (Gaussian) broadening in Domega [5.].")
   parser.add_option("-t", "--test", action="store_true",
                     help="Run internal doctests")
   parser.add_option("-n", "--dont_plot", action="store_true",
                     help="Do not plot results")
   parser.add_option("-p", "--plot", action="store_true",
                     help="Directly plot results (default if no output)")
   parser.add_option("-f", "--fast", action="store_true",
                     help="Tweak array sizes to low prime factors for fft")
   parser.add_option("--ignore", type=float, default=0.,
                     help="Ignore initial IGNORE ps")

   options, args = parser.parse_args()

   if options.test:
      testmod()
      sys.exit(0)

   ac_name, ft_name, grep_name = None, None, None
   if len(args) > 4:
      parser.error("Need at most three positional argument")
   if len(args) >= 4:
      grep_name = args[3]
   if len(args) >= 3:
      ac_name = args[2]
   if len(args) >= 2:
      ft_name = args[1]
   if len(args) >= 1:
      in_ = open(args[0], "r")
   else:
      in_ = sys.stdin

   try:
      if options.dont_plot or (len(args) > 1 and not options.plot):
         raise ImportError
      import warnings
      warnings.filterwarnings('ignore', category=DeprecationWarning)
      import matplotlib.pyplot as plt
      warnings.filterwarnings('default', category=DeprecationWarning)
      fig = plt.figure()
      ax1 = fig.add_subplot(2,1,1)
      ax2 = fig.add_subplot(2,1,2)
   except ImportError:
      plt = None
      ax1 = None
      ax2 = None


   sys.stderr.write("Reading input...\n")
   F, Dt = read_aims(in_, options.velocities, grep_name)
   if in_ is not sys.stdin:
      in_.close()

   n_timesteps, ndf = F.shape
   if n_timesteps == 0:
      error("Quantity not found")
   else:
      sys.stderr.write("Found %i time steps of %s fs. -> T_tot = %s ps\n" %
                       (n_timesteps, Dt*1000., n_timesteps*Dt))


   if options.ignore > 0.:
      n_ignore = int(options.ignore / Dt)
      sys.stderr.write("Ignoring the first %i time steps\n" % (n_ignore))
      F = F[n_ignore:, :]
      n_timesteps -= n_ignore

   # === Remove mean values

   # There may be a finite mean dipole moment, which can lead to a huge
   # zero-frequency mode.  Remove this one.
   for i in xrange(ndf):
      mean = np.mean(F[:,i])
      F[:,i] -= mean


   n_old = n_timesteps
   if options.fast:
      n_timesteps = to_high_two(n_timesteps)
      sys.stderr.write("Using only %i time steps to enable FFT."
                       "-> T_use = %s ps\n" %
                       (n_timesteps, n_timesteps*Dt))
      F = F[-n_timesteps:, :]

   # === Calculate autocorrelation

   N = int((1.-options.cutoff_ratio) * n_timesteps)
   N_old = N
   if options.fast:
      N = to_high_two(N)

   sys.stderr.write("Getting autocorralation from %s ps to %s ps...\n" %
                    (-N*Dt, N*Dt))
   f = autocorr(F, N)
   f_cut = f * Hanning(N)
   f_cut_ext = extend_to_whole(f_cut)
   
   # === Multiply with Hanning window


   if ac_name is not None:
      ac_out = open(ac_name, "w")
      ac_out.write("# Autocorrelation function f(t)\n")
      ac_out.write("# %10s %12s %12s\n" % (
            "t (ps)", "f(t)", "f(t)*Hanning(t)"))
      for i in xrange(N):
         ac_out.write("%12g %12g %12g\n" % (i*Dt, f[i], f_cut[i]))
      ac_out.close()

   if ax1 is not None:
      t = np.arange(-N+1, N) * Dt
      fe = extend_to_sym(f)
      fce = extend_to_sym(f_cut)
      ax1.plot(t, fe, label="Autocorrelation")
      ax1.plot(t, fce, label="Autocorrelation x Hanning")
      ax1.set_xlabel(r'$\Delta t$ (ps)')

   # === Perform fft

   sys.stderr.write("Performing Fourier transform...\n")
   ft_cmplx = np.fft.rfft(f_cut_ext)
   # Imaginary part should vanish for a symmetric original array:
   if not np.allclose(np.imag(ft_cmplx), 0.):
      sys.stderr.write("Maximum imaginary part in fft: %g\n" %
                       np.max(np.imag(ft_cmplx)))
   ft = np.real(ft_cmplx)
   assert len(ft) == N+1 # ft[N] contains Nyquist frequency.
   
   # === Broaden

   # broadening does not work yet
   T = (2*N-1) * Dt
   Domega = 2.*np.pi / T
   Dwave = Domega / C_LIGHT_CM_PS / (2*np.pi)
   sys.stderr.write("Correlation time of %.2f ps leads to"
                    " frequency steps of %.2f ps^-1 or %.2f cm^-1.\n"
                    % (T, Domega, Dwave))
   sigma = options.broadening # / Dwave
   n_gauss = min(N, 5 + 10*int(sigma))
   n_gauss = 2*(n_gauss/2) - 1 # make it odd
   gft = np.convolve(ft, Gaussian(n_gauss, sigma), 'same')

   # === Output

   if ft_name is not None:
      ft_out = open(ft_name, "w")
      ft_out.write("# Fourier transform of autocorrelation function\n")
      ft_out.write("# nu (cm^-1)\t f(nu)\t G*f(nu)\n")
      for i in range(N):
         ft_out.write("%g\t%g\t%g\n" % (i*Dwave, ft[i], gft[i]))
      ft_out.close()

   if ax2 is not None:
      waves = np.arange(N+1) * Dwave
      f_ext = extend_to_whole(f[:N])
      t = extend_to_whole(np.arange(N))
      ax2.plot(waves, ft, label="Fourier transform", color='#A0A0FF')
      ax2.plot(waves, gft, label="Broadend transform")
      ax2.set_xlabel(r'$\nu$ (cm$^{-1}$)')

      max_gft = np.max(gft)
      thres = 0.01 * max_gft
      for i in range(N-1, 0, -1):
         if gft[i] > thres: break
      n = min(i + 20, N)
      ax2.set_xlim([0, n*Dwave])


   if plt is not None:
      def on_q_exit(event):
          if event.key == "q": sys.exit(0)
      plt.figure(1).canvas.mpl_connect('key_press_event', on_q_exit)
      plt.show()

if __name__ == "__main__":
   _main()
