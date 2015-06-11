from math import *
import numpy as np
import matplotlib.pyplot as plt
from Functions import MasterEquationSolver_6levels
import sys

# DESCRIPTION
# Solve master equation for a system (6-levels) + vibronic-system coupling

# Five eigenbasis states are arranged in the folowing way :
# |0> = |d> == donor
# |1> = |b> ==  bridge 1
# |2> = |b> == bridge 2
# |3> = |a> == acceptor (para position)
# |4> = |b> == bridge 3
# |5> = |b> == bridge 4
# |R> for tunneling state from acceptor cite

# Hamiltonian
# H = \omega_0/2 ( |d><d| + |a><a| + \Sum_i |b_i><b_i| ) +
# \beta \Sum_i (|i><i+1| +hc)





# NOTE: 1eV = 2.41798e14 Hz

def main(argv):
    
    # Plot number
    plnum = argv[0]
    
    # Save plots true or false
    save = 1

    
    print "*************************************************"
    print "Coherence property of benzene"
    print "*************************************************"


    # BEGIN PROGRAM
    num = 6  # total number of basis states (without |R> state)

    # Acceptor position
    # 1 == ortho, 2 == meta, 3 == para (python numbering starts at 0)
    acptr_rowind = 1
    acptr_colind = acptr_rowind
    print "acceptor position : %d" % acptr_rowind


    if acptr_rowind==1:
        name =  "ortho"
    elif acptr_rowind==2:
        name = "meta"
    elif acptr_rowind==3:
        name = "para"






    # Constants

    # Planck's Constants (J.s)
    h = 6.626068e-34

    # Reduced Planck's Constant
    hbar = h/2/pi

    # Boltzmann's Constant (Joules/Kelvin)
    kb = 1.3806503e-23

    # NOTE: 1eV = 2.41798e14 Hz

    # System parameters

    w0 = 2*pi*1e13 # cite energy alpha (donor, acceptor, bridge are all equal)
    alpha = w0
    print "omega (or alpha) : %d" % w0


    beta = w0  # hopping
    print "beta : %d" % beta


    # Electron decay rate at acceptor cite
    Gamma_ed = w0*100.0
    print "Gamma_ed : %d" % Gamma_ed



    # Phonon bath parameters

    # Damping rate (in Hz) (phenomenological)
    Gamma = 1.0/10.0*w0
    print "Gamma : %d" % Gamma


    # Temperature (in Kelvin)
    Tempr = 300.0   # 290 to 300 K == room temperature
    print "Temperature : %d"  % Tempr


    # Time

    # True time (in s)
    runtime = 2*1e-12

    # Time Span (unitless parameter or in the units of 1/w0)
    tfinal = runtime*w0
    timeSpan = np.arange(1e-1, tfinal, 0.5) # 1D array
    timeSpan = np.reshape(timeSpan, (timeSpan.size, )) # 1D array


    # Calculate time dependent density matrix rho(t) numerically
    # Call MasterEquationSolver

    # Input format of rho matrix
    # Vertically concantenating rows of a density matrix
    # Thus rho0 = [rho11; rho12; rho13; rho14;.....; rho21; rho22;
    # rho23; rho24; rho24;.......; rho41; rho42; rho43; rho44;...]

    # NOTE : rho(0) = column vector

    # Output format of rho matrix
    # rho(t) = [rho11, rho12, rho13, rho14,....., rho21, rho22,
    # rho23, rho24, rho24,......., rho41, rho42, rho43, rho44,...]

    # NOTE : rho(t) = matrix
    # matrix whose each row is values of rho at different time t (Note: the
    # first row is rho(0) column vector)

    # Intial condition
    rho0Mat= np.zeros((num, num), dtype=np.complex128)
    rho0Mat[0, 0] = 1.0

    # Reshape initial density matrix to a column vector
    rho0real = np.concatenate((np.reshape(np.real(rho0Mat), (num**2, )), np.zeros((1, ))), axis=0)  # with state |R>
    rho0img = np.concatenate((np.reshape(np.imag(rho0Mat), (num**2, )), np.zeros((1, ))), axis=0)  # with state |R>
    rho0 = np.concatenate((rho0real, rho0img), axis=0) # with state |R>

    # rho0real = np.reshape(np.real(rho0Mat), (num**2, )) # without state |R>
    # rho0img = np.reshape(np.imag(rho0Mat), (num**2, )) # without state |R>
    # rho0 = np.concatenate((rho0real, rho0img), axis=0) # without state |R>
    # 1D array


    # Calculate rho(t)
    [rho, timeSpan] = MasterEquationSolver_6levels(timeSpan, rho0, kb,  hbar, w0, alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed, Tempr, num)


    if rho0.size==2*(num**2+1): # with |R> state
        rho = rho[:, 0:num**2+1] + 1j*rho[:, num**2+1:]
    elif rho0.size==2*(num**2): # without |R> state
        rho = rho[:, 0:num**2] + 1j*rho[:, num**2:]




    # PLOTS

    plotind = np.diag(np.reshape(np.arange(0, num**2), (num, num)), 0)  # diagonal entires (read only)
    nrows = plotind.size   # Counts number of indices in plotind

    # Make Strings
    rholabels = []
    for i in range(0, nrows, 1):
        (m, n) = np.unravel_index(plotind[i], (num, num)) # matlab ind2sub
        # Python row major vs Matlab column major
        rholabels.append(r"$\rho_{" + str(m) + str(n) + "}$" + "(t)")


    # Finally plot them

    # Number of plots per figure (ppf)
    ppf = 3

    quotient = nrows/ppf
    for n in range(0, quotient, 1):
        fn = plt.figure()
        for i in range(0, ppf, 1):
            # Note : subplot(nrows, ncolumn, plot_numb>0)
            plt.subplot(ppf, 1, i+1)
            plt.plot(timeSpan/w0*1e12,  np.real(rho[:, plotind[i+n*ppf]]), "r")
            plt.ylabel(rholabels[i+n*ppf])
            plt.draw()
        plt.xlabel("Time (ps)")
        # Save
        if save:
            plt.savefig("plots6/" + "rho_a" + str(n) + "_" + str(plnum), format="pdf")



    rmdr = nrows % ppf  # or np.remainder(nrows, ppf)
    if rmdr!=0:
        fr = plt.figure()
        for i in range(0, rmdr, 1):
            # Note : subplot(nrows, ncolumn, plot_numb)
            plt.subplot(rmdr, 1, i)
            plt.plot(timeSpan/w0*1e12,  np.real(rho[:, plotind[i+quotient*ppf]]), "r")
            plt.ylabel(rholabels[i+quotient*ppf])
            plt.draw()
        plt.xlabel("Time (ps)")
        # Save
        if save:
            plt.savefig("plots6/" + "rho_a" + str(quotient) + "_" + str(plnum), format="pdf")





    # Plot survival probability

    fs = plt.figure()
    ax = fs.add_subplot(1,1,1)
    plt.plot(timeSpan/w0*1e12, np.sum(np.real(rho[:, plotind]), axis=1)) # sum returns 1D array
    plt.xlabel("Time (ps)")
    plt.ylabel("Survival probability")
    plt.text(0.2, 0.9, name, fontsize=12, horizontalalignment="left", verticalalignment="top",
             transform=ax.transAxes)
    plt.draw()


    # Save
    if save:
        plt.savefig("plots6/" + "Psurv_" + str(plnum), format="pdf")


    # Show plots
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:]) #sys.argv[0] is the name of the program file
