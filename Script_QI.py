

# DESCRIPTION

# Full Hamiltonian of the system is
# H = H_s + H_ph + H _{el-ph}
# Here we will assume that each electronic system is coupled to individual
# phonon bath


# Electronic system Hamiltonian in the second-quantized notation:
# H_s = \alpha \sum_i c^\dag_i c_i + \beta \sum_i (c^\dag_i c_{i+1} + h. c.)

# The phonon baths Hamiltonian is given by:
# H_ph = \sum_i \sum_k \hbar \omega_{i, k} b^\dag_{i, k} b_{i, -k}


# The electron-phonon coupling Hamiltonian is given by:
# H_{el-ph} \sum_i \sum_k \hbar \g_{i, k} [c_i b^\dag_{i, k} + c^\dag_i b_{i, -k}]
# + \sum_i \sum_k \hbar \f_{i, k} c\dag_i c_i[ b^\dag_{i, k} + b_{i, -k}]



# We focus on zero excitation (or ground state) and single
# electron/excitation subspace of full Hilbert space:

# In zero and single excitation subspace, the basis sets are arranged in
# the following way:

# |L> = left lead connected to the donor (site 1)
# |G> = ground state with zero excitation
# |1> = single excitation at site 1
# |2> = single excitation at site 2
# |3>  = ...
# |4> 
# |5> 
# |6> = single excitation at site 6
# |R>  = right lead connected to the acceptor (any one of the sites other than site 1)

# NOTE: 1eV  <--> 2.41798e14 Hz (linear frequency)



 
# List all available backends (1st step)
# import matplotlib
# print matplotlib.rcsetup.all_backends
 
import matlab.engine  # This must be imported before importing user defined functions/objects
eng = matlab.engine.start_matlab()

 
# Choose backend (2nd step)
import matplotlib
matplotlib.use('TkAgg') # matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt    
    
# Get the backend in use (3rd step)
# print matplotlib.pyplot.get_backend()
    
    
import sysconfig, sys, os, time

from math import *
import numpy as np
from scipy import linalg


from joblib import Parallel, delayed
import multiprocessing

from multiprocessing import Pool

# Import user defined functions/objects
from Functions import MasterEquationSolver_QI, PlankFun


    

    
    
# Converts 2D numpy array type to Matlab matrix type
def N2M(A):
    
    A = np.array(A)
    
    tup = np.shape(A)
    
    if len(tup)<2:
        r = tup[0]
        c = 1
    else:
        r = tup[0]
        c = tup[1]
    
    # flatten numpy 2D array to a list
    B = A.flatten('C').tolist() # row major or C-style
    
    # Converts to matlab type
    Bm = matlab.double(B)
    Bm = eng.reshape(Bm, r, c)
    
    return Bm
        

def main(argv):
    
    
    # Start time stamp
    date_string = time.asctime() 
    print '\n  START TIME: % s \n' % date_string 

    # Python or Matlab for plotting
    Python = 0
    Matlab = 1
        
   
    # Number of cores
    num_cores = os.environ.get('NUMCORES')  # multiprocessing.cpu_count() # 
    if num_cores==None:
        num_cores = 1 
    else:
        num_cores = int(num_cores)
        
    print 'number of cores: %d' %num_cores
    
    
    # Name of the folder
    folder = os.environ.get('FOLDERNAME') # 'plots' 
    if not os.path.isdir(folder):
        os.mkdir(folder, 0755)
            
    # Plot number
    plnum = int(argv[0])


    # Run calculation (true or false)
    runcalc = int(argv[1])
    
    # Save plots true or false
    savefig = int(argv[2])
    
    
    print "********************************"
    print "Transport property of benzene"
    print "*********************************"
    
    

    # BEGIN PROGRAM
    num_sites = 6
    num = num_sites + 3  # total number of basis states 



    # Constants
    h = 6.626068e-34 # Planck's Constants (J.s)
    eps0 = 8.812*1e-12 
    c = 3*1e8 
    hbar = 1.05*1e-34 
    me = 9.109*1e-31 
    kb = 1.3806503e-23 #  Boltzmann's Constant (Joules/Kelvin)
    ec = 1.6e-19 
    aB = 5.3e-11
    ev2hz = ec/2.0/pi/hbar
    hz2ev = 1.0/ev2hz
    invcm2hz = 100.0*c
    hz2invcm = 1.0/invcm2hz


    # Donor, left and right leads numbering 
    # 0 = L
    # 1 = G
    # 2 == donor

    # Acceptor numbering
    # 3 == ortho, 4 == meta, 5 == para (python numbering starts from 0)
    acptr_ind = 3
    print "acceptor position : %d" % acptr_ind


    if acptr_ind==3:
        name =  "ortho"
    elif acptr_ind==4:
        name = "meta"
    elif acptr_ind==5:
        name = "para"



    # System parameters

    alpha = 2*pi*1e13 # site energy alpha (donor, acceptor, bridge are all equal)
    w0 = alpha 
    print "alpha (or w0) : %d" % alpha


    beta = w0  # hopping
    print "beta : %d" % beta



    # Phonon bath parameters

    # Damping rate (in rad/sec) (phenomenological)
    Gamma = 0.0/10.0*w0
    print "Gamma : %d" % Gamma


    # Dephasing rate 9in rad/sec)
    Gamma_phi = 0.0/10.0*w0
    print "Gamma_phi : %d" % Gamma_phi


    # Temperature (in Kelvin)
    Tempr = 300.0   # 290 to 300 K == room temperature
    print "Temperature : %d"  % Tempr



    # Electron hopping from left and right leads (in rad/sec)
    
    Gamma_donor = w0*100.0 
    print "Gamma donor : %d" % Gamma_donor


    Gamma_acptr = w0*100.0
    print "Gamma acceptor : %d" % Gamma_acptr



    # System Hamiltonian

    # Linear index of acceptor site
    acptr_linind = np.ravel_multi_index((acptr_ind, acptr_ind), (num, num))
    # matlab sub2ind (Python has row major. Matlab has column major)
    

    # System Hamiltonian (in site basis {|i>})
    H_site = alpha*np.diag(np.ones(num_sites), 0)
    H_hop = beta*np.diag(np.ones(num_sites-1), -1) + beta*np.diag(np.ones(num_sites-1), 1)
    H_hop[0, num_sites-1] = beta
    H_hop[num_sites-1, 0] = beta
    # benzene ring (periodic boundary)

    # Take left lead |L>, ground state |G> and right lead |R> states into account
    H_site = linalg.block_diag(linalg.block_diag(np.zeros((2, 2)), H_site), np.zeros((1, 1)))
    H_hop =  linalg.block_diag(linalg.block_diag(np.zeros((2, 2)), H_hop), np.zeros((1, 1)))

    # Total Hamiltonian
    H = hbar/2.0*(H_site + H_hop)


    

    # Jump operators

    # Collapse operators for each site (scaled by square root of decay rate)
    ScaledCollapseOperators = []
    A = np.zeros((num, num), dtype=np.float64)
    for i in range(2, num-1):
        A[1, i] = 1.0
        Adag = np.conjugate(A).T
        ScaledCollapseOperators.append(np.sqrt(Gamma*(1.0+PlankFun(alpha, Tempr)))*A)
        ScaledCollapseOperators.append(np.sqrt(Gamma*PlankFun(alpha, Tempr))*Adag)
        ScaledCollapseOperators.append(np.sqrt(Gamma_phi)*np.dot(Adag, A))

    # Collapse operator for input excitation to donor
    # (scaled by square root of coupling rate of donor to the left lead)
    B = np.zeros((num, num), dtype=np.float64)
    B[0, 2] = 1.0
    Bdag = np.conjugate(B).T
    ScaledCollapseOperators.append(np.sqrt(Gamma_donor)*Bdag)

    # Collapse operator for output excitation from acceptor 
    # (scaled by square root of coupling rate of acceptor to the right lead)
    C = np.zeros((num, num), dtype=np.float64)
    C[-1, acptr_ind] = 1.0
    ScaledCollapseOperators.append(np.sqrt(Gamma_acptr)*C)




    # Time

    # True time (in s)
    runtime = 20.0/w0

    # Time Span (unitless parameter or in the units of 1/w0)
    tfinal = runtime*w0
    tspan = np.linspace(1e-3, tfinal, 1000) # 1D array (not a column vector of a 2D array)
    tspan = np.reshape(tspan, (tspan.size, )) # 1D array (not a column vector of a 2D array)


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
    rho0Mat[0, 0] = 1.0 # at t = 0, left lead is occupied
    

    # Reshape initial density matrix to a column vector
    rho0 = np.reshape(rho0Mat, (num**2, ))
    

    # Calculate rho(t)
    [rho, timespan] = MasterEquationSolver_QI(tspan, rho0, H, ScaledCollapseOperators, w0, num)
    
    



    # **************************
    
    # PLOTS
    if savefig==1:

        
        print 'begin plotting ... '
        indx = np.reshape(np.arange(0, num**2), (num, num)) # C-major or row major
        plotind = np.diag(indx, 0)  # diagonal entries (read only)            
        # plotind = np.concatenate((plotind, indx[np.triu_indices(num, 1)]), axis=0) # off-diagonal entries
        nrows = plotind.size   # Counts number of indices in plotind

        # Make Strings
        rholabels = []
        for i in range(0, nrows, 1):
            (m, n) = np.unravel_index(plotind[i], (num, num)) # matlab ind2sub
            # Python row major vs Matlab column major
            rholabels.append("$\rho_{" + str(m) + str(n) + "}$" + "(t)")


        # Finally plot them
    

        # Number of plots per figure (ppf)
        ppf = 3

        quotient = nrows/ppf
        for n in range(0, quotient, 1):
            folder_filename = folder + "/rho_a" + str(n) + "_" + str(plnum)        
            if Python==1:
                fn = plt.figure()
            if Matlab==1:    
                fn = eng.figure(nargout=0)
            for i in range(0, ppf, 1):
                if Python==1:
                    # Note : subplot(nrows, ncolumn, plot_numb>0)
                    plt.subplot(ppf, 1, i+1)
                    plt.plot(tspan/w0*1e12,  np.abs(rho[:, plotind[i+n*ppf]]), "r")
                    #plt.ylabel(rholabels[i+n*ppf])
                    plt.draw()
                if Matlab==1:
                    eng.subplot(float(ppf), 1.0, i+1.0, nargout=0)
                    eng.plot(N2M(tspan/w0*1e12),  N2M(np.abs(rho[:, plotind[i+n*ppf]])), 'r')
                    eng.ylabel(rholabels[i+n*ppf])
                    eng.hold('off', nargout=0)
            if Python==1:
                plt.xlabel("Time (ps)")
                plt.savefig(folder_filename, format="pdf")
            if Matlab==1:
                eng.xlabel('Time (ps)')
                eng.savefig(folder_filename + '.fig', nargout=0)
                    

        rmdr = nrows % ppf  # or np.remainder(nrows, ppf)
        if rmdr!=0:
            if Python==1:
                fr = plt.figure()
            if Matlab==1:
                fr = eng.figure(nargout=0)
            for i in range(0, rmdr, 1):
                if Python==1:
                    # Note : subplot(nrows, ncolumn, plot_numb)
                    plt.subplot(rmdr, 1, i+1)
                    plt.plot(tspan/w0*1e12,  np.abs(rho[:, plotind[i+quotient*ppf]]), "r")
                    # plt.ylabel(rholabels[i+quotient*ppf])
                    plt.draw()
                if Matlab==1:
                    eng.subplot(float(rmdr), 1.0, i+1.0, nargout=0)
                    eng.plot(N2M(tspan/w0*1e12),  N2M(np.abs(rho[:, plotind[i+quotient*ppf]])), 'r')
                    eng.ylabel(rholabels[i+quotient*ppf])
                    eng.hold('off', nargout=0)
            
            folder_filename = folder + "/rho_a" + str(quotient) + "_" + str(plnum)
            if Python==1:
                plt.xlabel("Time (ps)")
                plt.savefig(folder_filename, format="pdf")
            if Matlab==1:
                eng.xlabel('Time (ps)')    
                eng.savefig(folder_filename + '.fig', nargout=0)

        
        # Plot survival probability
        folder_filename = folder + "/P_surv_escp_" + str(plnum)
        if Python==1:
            fp = plt.figure()
            plt.plot(tspan/w0*1e12,  np.sum(np.abs(rho[:, plotind[2:-1]]), 1), "r")
            plt.plot(tspan/w0*1e12,  np.abs(rho[:, plotind[-1]]), "b")
            plt.draw()
            plt.xlabel("Time (ps)")
            # plt.ylabel(rholabels[i+quotient*ppf])
            plt.savefig(folder_filename, format="pdf")
        if Matlab==1:
            fp = eng.figure(nargout=0)
            eng.plot(N2M(tspan/w0*1e12),  N2M(np.sum(np.abs(rho[:, plotind[2:-1]]), 1)), "r")
            eng.hold('on', nargout=0)
            eng.plot(N2M(tspan/w0*1e12),  N2M(np.abs(rho[:, plotind[-1]])), "b")
            eng.hold('off', nargout=0)
            eng.xlabel("Time (ps)")
            eng.ylabel("Probability")
            eng.savefig(folder_filename + '.fig', nargout=0)

    # Show plots
    if Python==1:
        plt.show()
        
    # Quit matlab    
    if Matlab==1:    
        print 'matlab stopping..'
        eng.quit()
    
    # End time stamp
    date_string = time.asctime() 
    print '\n  END TIME: % s \n' % date_string





if __name__ == "__main__":
    # Home
    os.environ['NUMCORES'] = '4'
    os.environ['FOLDERNAME'] = 'plots_QI'
    argv = [1, 1, 1]
    main(argv)
    
    # Cluster
    # main(sys.argv[1:]) #sys.argv[0] is the name of the program file





# Put text in  a plot
# plt.text(0.2, 0.9, name, fontsize=12, horizontalalignment="left", verticalalignment="top",\
             # transform=ax.transAxes)
