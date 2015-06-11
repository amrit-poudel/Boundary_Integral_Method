import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scpintg
from math import *




# Derivative of rho
def rhoPrime_6levels(t, rho, kb,  hbar, w0, alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed, Tempr, num):

    # Description
    # Evaluates derivatives of rho vector

    # Master equation takes the form: drho/dt =  D/w0 , where t = untiless (units of 1/w0)
    # D = derivative of "rho" vector and "rho" is a
    # vector formed from the density matrix "rho"  in the following way:
    # Vertically concantenating rows of a density matrix
    # Thus rho = [rho11; rho12; rho13; rho14;.....; rho21; rho22;
    # rho23; rho24; rho24;.......; rho41; rho42; rho43; rho44;...]

    # Input
    # t = time
    # rho = column vector of elements of density matrix
    # kb, hbar, w0, alpha, beta, Gamma, Gamma_ed, Tempr are several parameters defined in Script
    # num = total number of basis states

    # Output
    # D = Derivative of a rho vector (D has units of frequency, that is, Hz)



    # Pauli operators
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    sy = np.array([[0.0, -1j],[1j, 0.0]])
    sz = np.array([[1.0, 0.],[0.0, -1.0]])
    sp = (sx + 1j*sy)/2.0
    sm = (sx - 1j*sy)/2.0

    # Linear index of acceptor cite
    acptr_linind = np.ravel_multi_index((acptr_rowind, acptr_colind), (num, num))
    # matlab sub2ind (Python has row major. Matlab has column major)
    

    # System Hamiltonian (in site basis {|i>})
    H_cite = w0*np.diag(np.ones(num), 0)
    H_hop = beta*np.diag(np.ones(num-1), -1) + beta*np.diag(np.ones(num-1), 1)
    H_hop[0, num-1] = beta
    H_hop[num-1, 0] = beta
    # benzene ring (periodic boundary)



    # Total Hamiltonian
    H = hbar/2*(H_cite + H_hop)

    (EE, EV) = np.linalg.eig(H)

    wc = np.max(EE)/hbar # cutoff frequency
    
    # Begin computing rhodot

    # Convert a rho column vector into rho matrix
    if rho.size==2*(num**2+1): # with state |R>
        rhoreal = rho[0:num**2]
        rhoimg = rho[num**2+1:-1]
        rhoMat = np.reshape(rhoreal, (num, num)) + 1j*np.reshape(rhoimg, (num, num))
    elif rho.size==2*num**2: # without state |R>
        rhoreal = rho[0:num**2]
        rhoimg = rho[num**2:]
        rhoMat = np.reshape(rhoreal, (num, num)) + 1j*np.reshape(rhoimg, (num, num))



    # This uses phenomenological model for dissipation
    rhodot =  1/1j/hbar*(np.dot(H, rhoMat) - np.dot(rhoMat, H))


    for m in range(0, num, 1):
        for n in range(0, num, 1):
            w = (EE[m] - EE[n])/hbar
            if np.abs(w)>1:
                # Lindblad operators (defined in energy eigenbasis and trasformed back to site basis {|i>} using energy eigenvectors)
                A = np.dot(EV[:, n], EV[m, :])
                Adag = np.conj(A).T  #  or np.conjugate(A).T or A.conjugate().T or A.conj().T
    
                rhodot = rhodot + \
                Gamma*pi*w*wc/(w**2 + wc**2)*(1.0 + 1.0/(exp(hbar*w/kb/Tempr) - 1.0))*(np.dot(np.dot(A, rhoMat), Adag) - \
                1/2.0*np.dot(np.dot(Adag, A), rhoMat) - 1/2.0*np.dot(np.dot(rhoMat, Adag), A))


#    # This approach is not quite correct since Lindblad operators are defined in site  basis {|i>} instead of energy eigenbasis
#    for m in range(0, num, 1):
#    
#        # Lindblad operators  (in site basis {|i>})
#        A = np.zeros((num, num))
#        
#        n = m + 1
#        if n>num-1:
#            n = 0
#        A[m, n] = 1.0
#        Adag = np.conj(A).T  #  or np.conjugate(A).T or A.conjugate().T or A.conj().T
#        
#        rhodot = rhodot + \
#        Gamma*(1.0 + 1.0/(exp(hbar*w0/kb/Tempr) - 1))*(np.dot(np.dot(A, rhoMat), Adag) - \
#        1/2.0*np.dot(np.dot(Adag, A), rhoMat) - 1/2.0*np.dot(np.dot(rhoMat, Adag), A)) + \
#        Gamma*1.0/(exp(hbar*w0/kb/Tempr)-1.0)*(np.dot(np.dot(Adag, rhoMat), A) - \
#        1/2.0*np.dot(np.dot(A, Adag), rhoMat) - 1/2.0*np.dot(np.dot(rhoMat, A), Adag))
#




    # Include tunneling to the right state  (R)
    rhodot[acptr_rowind, acptr_colind] = rhodot[acptr_rowind, acptr_colind] - Gamma_ed*rhoMat[acptr_rowind, acptr_colind]  # diagonal
            
    for i in range(0, num, 1):
        if (i!=acptr_colind):
            rhodot[i, acptr_colind] = rhodot[i, acptr_colind] - Gamma_ed/2*rhoMat[i, acptr_colind]   # off-diagonal
        
        if (i!=acptr_rowind):
            rhodot[acptr_rowind, i] = rhodot[acptr_rowind, i] - Gamma_ed/2*rhoMat[acptr_rowind, i]   # off-diagonal


    # Output rhodot in a vector form
    D = np.reshape(rhodot, (num**2, ))  # 1D array
    
    
    # Include derivative of the right  state due to tunneling from the acceptor cite
    if rho.size==2*(num**2+1):
        D = np.concatenate((D, np.ones(1, )*Gamma_ed*rhoMat[acptr_rowind, acptr_colind]), axis=0) # 1D array




    return np.concatenate((np.real(D)/w0, np.imag(D)/w0), axis=0)




# Master equation solver
def MasterEquationSolver_6levels(tspan, rho0, kb,  hbar, w0, alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed, Tempr, num):
    
    # Description
    # Solves Master Equation represented in energy eigenfunction basis
    # Master equation takes the form: drho/dt =  D/w0 , where t = untiless (units of 1/w0) and "rho" is a
    # vector formed from the density matrix "rho"  in the following way:
    # rho0 = [rho11; rho12; rho13; rho14;.....; rho21; rho22;
    # rho23; rho24; rho24;.......; rho41; rho42; rho43; rho44;...]
    
    # NOTE : rho(t) = matrix
    
    # Input
    # tpan = time span
    # rho0 = rho(0), represented as a column vector
    # kb, hbar, w0,alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed, Tempr  are several parameters defined in Script
    # num = total number of basis states
    
    # Output
    # rho = matrix whose each row is values of rho at different time t (Note: the
    # first row is rho(0) column vector )
    
    reltol = 1e-8*np.ones(rho0.shape)
    abstol = 1e-7*np.ones(rho0.shape)
    
    # rho = scpintg.odeint(rhoPrime_6levels, rho0, tspan, args=(kb, hbar, w0, alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed,
    # Tempr, num), rtol=reltol, atol=abstol)
                                                          
    solver = scpintg.ode(rhoPrime_6levels).set_integrator('vode', method='bdf', atol=abstol, rtol=reltol)
    solver.set_initial_value(rho0, tspan[0])
    solver.set_f_params(kb, hbar, w0, alpha, beta, acptr_rowind, acptr_colind, Gamma, Gamma_ed,
    Tempr, num)
    

    num_steps = np.size(tspan)
    delta_t = tspan[1]-tspan[0]

    rho = np.zeros((num_steps, np.size(rho0)))
    rho[0, :] = rho0
    
    t = np.zeros((num_steps,))
    t[0] = tspan[0]
    
    k = 1
    while k < num_steps: # solver.successful() and
        solver.integrate(solver.t + delta_t)
        
        # Store the results
        t[k] = solver.t
        rho[k, :] = solver.y
        k += 1

    return rho, t


