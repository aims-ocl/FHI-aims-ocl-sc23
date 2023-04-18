#!/usr/bin/python

"""
The parameters in the Birch-Murnaghan equation of state

  E(V) = E0 + B0*V/B' * [(V0/V)^B' / (B'-1) + 1] - B0 V0 / (B'-1)

are fit to energies for different lattice constants.

  E0: Energy in equilibrium
  V0: Volume in equilibrium
  B0: Bulk modulus in equilibrium
      B := - V (dP/dV)_T
      P := - (dE/dV)_S
  B': Pressure derivative of B  (assumed constant in this model)
      B' := (dB/dP)_T

See also: http://en.wikipedia.org/wiki/Birch-Murnaghan_equation_of_state
"""

import sys

import numpy as np
from scipy.optimize import fmin_powell

USAGE = """%prog [options] <input_file>

The <input_file> should contain two columns, the first containing the volume
and the second total energies.  

No unit conversions are done. If you have used Angstroms^3 and eV in the input
file, the bulk modulus is returned in eV/Anstrom^3.  You can convert it using
GNU units:

  units "<B0> eV/angstrom^3" "GPa"\
"""

def Birch_Murnaghan(paras, V):
    """Evaluate Birch Murnaghan equation-of-state for [E0, V0, B0, B']"""
    E0, V0, B0, Bprime = paras
    T1a = B0*V/Bprime
    T1b = ((V0/V)**Bprime / (Bprime-1.) + 1.)
    T1 = T1a * T1b
    T2 = - B0 * V0 / (Bprime - 1)
    E = E0 + T1 + T2
    return E

def sum_of_squares(paras, Vs, Es):
    """Calculate squared norm of residue vector for paras."""
    E_BMs = Birch_Murnaghan(paras, Vs)
    return np.sum((Es - E_BMs)**2)


def initial_guess(Vs, Es):
    """Construct initial guess for V0, E0, ..."""
    n = len(Es)
    assert n == len(Vs)
    if n < 3:
        parser.error("Need at least three energies for fit")  # 4?

    E0 = np.min(Es)
    i0 = np.argmin(Es)   # index of minimum
    V0 = Vs[i0]
    if i0 > 0:
        dE_left = (Es[i0] - Es[i0-1]) / (Vs[i0] - Vs[i0-1])
        V_left = (Vs[i0] + Vs[i0-1]) / 2.
    else:
        dE_left = 0.
        V_left = Vs[i0]
    if i0 < n-1:
        dE_right = (Es[i0+1] - Es[i0]) / (Vs[i0+1] - Vs[i0])
        V_right = (Vs[i0+1] + Vs[i0]) / 2.
    else:
        dE_right = 0.
        V_right = Vs[i0]
    ddE = (dE_left - dE_right) / (V_left - V_right)
    B0 = ddE * V0
    Bprime = 3.5  # Typical value according to [1]
    return np.array([E0, V0, B0, Bprime])


def get_default_VV(V0, Vs):
    """Get default set of lattice constants aa for plotting"""
    Vmin = min(np.min(Vs), V0)
    Vmax = max(np.max(Vs), V0)
    Vrange = Vmax - Vmin
    VV = np.linspace(Vmin - 0.15*Vrange, Vmax + 0.15*Vrange, 200)
    return VV

def plot(VV, EE, Vs, Es, V_text="$V$ in [a^3]"):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Plot input data points
    ax.plot(Vs, Es, 'o', label="Calculated energies")

    # Plot Birch-Murnaghan equation of states
    ax.plot(VV, EE, '-', label="Birch-Murnaghan fit")

    # Decorate and plot
    ax.legend()
    ax.set_xlabel(V_text)
    ax.set_ylabel("$E$ in [E]")
    def on_q_exit(event):
        if event.key == "q": sys.exit(0)
    plt.connect('key_press_event', on_q_exit)
    plt.show()


def main():
    import optparse

    # Parse command line
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-p", "--plot", action="store_true",
                      help="Plot E(a) curve using matplotlib")
    parser.add_option("-i", "--info", action="store_true",
                      help="Output info about Birch-Murnaghan e.o.s. & quit")
    parser.add_option("-r", "--range", action="store", nargs=3, type="float",
                      help="Output range <start> <stop> <nstep>")
    parser.add_option("-o", "--output", "--datafile", action="store",
                      help="Output fitted curve to OUTPUT")
    parser.add_option("-e", "--reference", action="store", type="float",
                      help="Reference energy", default=0.)
    parser.add_option("-l", "--lattice-constant", metavar="FACTOR",
                      help="""\
Assume <input_file> contains lattice constants a instead of volumes V
and calculate V=FACTOR*a^3. Use, e.g., FACTOR=0.25 for fcc.""")

    options, args = parser.parse_args()
    if options.info:
        print __doc__
        sys.exit(0)
    if len(args) != 1:
        parser.error("Need exactly one argument")
    input_file_name = args[0]
    
    use_lattconst = options.lattice_constant is not None
    if use_lattconst:
        unit_cell_frac = float(options.lattice_constant)

    # Read input data
    data = np.loadtxt(input_file_name)   # Input data as array
    if use_lattconst:
        lattconsts, Es = data.T              # Separate 1st and 2nd col
        Vs = unit_cell_frac * lattconsts**3
    else:
        Vs, Es = data.T              # Separate 1st and 2nd col
    Es -= options.reference

    # Sort input data by V
    VEs = sorted([zip(Vs, Es)])          # List of (V, E) pairs sorted by V
    Vs, Es = np.array(VEs).T             # Convert back to numpy arrays

    # Construct initial guesses
    paras = initial_guess(Vs, Es)

    # Optimize parameters
    paras, fopt, direc, n_iter, n_funcalls, warnflag \
        = fmin_powell(sum_of_squares, paras, (Vs, Es), full_output=True)
    E0, V0, B0, Bprime = paras

    # Final output
    print "Final residue: %g [E]^2." % fopt
    print "=== Final parameters (output units match input units)"
    outlist = []
    outlist.append(" E0  = %18.12g [E]" % E0)
    outlist.append(" V0  = %18.12g [a^3]" % V0)
    if use_lattconst:
        a0 = (V0 / unit_cell_frac)**(1./3.)
        outlist.append(" a0  = %18.12g [a]" % a0)
    outlist.append(" B0  = %18.12g [E/a^3]" % B0)
    outlist.append(" B0' = %18.12g" % Bprime)
    for line in outlist:
        print line

    if options.range is None:
        VV = get_default_VV(V0, Vs)
    else:
        start, stop, nstep = options.range
        VV = np.linspace(start, stop, nstep)

    EE = Birch_Murnaghan(paras, VV)


    # Optionally plot
    if options.output is not None:
        f = open(options.output, "w")
        for line in outlist:
            f.write("# %s\n" % line)
        if use_lattconst:
            np.savetxt(f, np.array([aa, VV, EE]).T)
        else:
            np.savetxt(f, np.array([VV, EE]).T)
        f.close()
    if options.plot:
        try:
            if use_lattconst:
                aa = (VV / unit_cell_frac)**(1./3.)
                plot(aa, EE, lattconsts, Es, "$a$ in [a]")
            else:
                plot(VV, EE, Vs, Es)
        except ImportError:
            parser.error("Cannot use --plot without matplotlib")

if __name__ == "__main__":
    main()
