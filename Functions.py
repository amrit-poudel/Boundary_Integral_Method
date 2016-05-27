import numpy as np
import scipy.integrate as scpintg
from math import *


import numpy as np
import scipy.integrate as scpintg
from math import *
from scipy import linalg
import gzip

# Derivative of rho
def MasterEquation(t, rho, H, ScaledCollapseOperators, w0, num):

    # Description
    # Evaluates derivatives of rho vector

    # Master equation takes the form: drho/dt =  D/w0 , where t = untiless (units of 1/w0)
    # D = derivative of "rho" vector and "rho" is a
    # vector formed from the density matrix "rho"  in the following way:
    # Vertically concantenating rows of a density matrix
    # Thus rho = [rho11; rho12; rho13; rho14;.....; rho21; rho22;
    # rho23; rho24; rho24;.......; rho41; rho42; rho43; rho44;...]

    # Input
    # t = dimensionless time
    # rho = column vector of elements of density matrix
    # H = system Hamiltonian (in Joules)
    # ScaledCollapseOperators = collapse operators A scaled by a factor of \sqrt{\Gamma*(1+N(w))}
    # w0 = characteristic frequency  in rad/sec(used to define time scale)
    # num = total number of basis states
    
    # Output
    # D = Derivative of a rho vector (D has units of frequency, that is, Hz)

    
    # Constants
    h = 6.626068e-34 # Planck's Constants (J.s)
    hbar = 1.05*1e-34 
    kb = 1.3806503e-23 #  Boltzmann's Constant (Joules/Kelvin)
    

    # Pauli operators
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    sy = np.array([[0.0, -1j],[1j, 0.0]])
    sz = np.array([[1.0, 0.],[0.0, -1.0]])
    sp = (sx + 1j*sy)/2.0
    sm = (sx - 1j*sy)/2.0
    

    (EE, EV) = np.linalg.eig(H)

    wc = np.max(EE)/hbar # cutoff frequency
    
    # Begin computing rhodot

    # Convert a rho column vector into rho matrix
    rhoreal = rho[0:num**2]
    rhoimg = rho[num**2:]
    rhoMat = np.reshape(rhoreal, (num, num)) + 1j*np.reshape(rhoimg, (num, num))



    # This uses phenomenological model for dissipation
    rhodot =  1.0/1j/hbar*(np.dot(H, rhoMat) - np.dot(rhoMat, H))

    
    # Lindblad equation
    for i, C in enumerate(ScaledCollapseOperators):
            rhodot = rhodot + LindbladFun(C, C, rhoMat)
    



    # Output rhodot in a vector form
    D = np.reshape(rhodot, (num**2, ))  # 1D array
    

    return np.concatenate((np.real(D)/w0, np.imag(D)/w0), axis=0)




# Master equation solver
def MasterEquationSolver_QI(tspan, rho0, H, ScaledCollapseOperators, w0, num):
    
    # Description
    # Solves Master Equation represented in energy eigenfunction basis
    # Master equation takes the form: drho/dt =  D/w0 , where t = untiless (units of 1/w0) and "rho" is a
    # vector formed from the density matrix "rho"  in the following way:
    # rho0 = [rho11; rho12; rho13; rho14;.....; rho21; rho22;
    # rho23; rho24; rho24;.......; rho41; rho42; rho43; rho44;...]
    
    # NOTE : rho(t) = matrix
    
    # Input
    # tpan = time span (dimensionless)
    # rho0 = column vector of rho(0)
    # H = system Hamiltonian (in Joules)
    # ScaledCollapseOperators = collapse operators A scaled by a factor of \sqrt{\Gamma*(1+N(w))}
    # w0 = characteristic frequency  in rad/sec(used to define time scale)
    # num = total number of basis states
    
    
    # Output
    # rho = matrix whose each row is values of rho at different time t (Note: the
    # first row is rho(0) column vector )
    
    # Convert complex rho0 2D array to a real 1D vector
    rho0real = np.reshape(np.real(rho0), (num**2, )) 
    rho0img = np.reshape(np.imag(rho0), (num**2, )) 
    rho_init = np.concatenate((rho0real, rho0img), axis=0) 
    # 1D array
    

    reltol = 1e-8*np.ones(rho_init.shape)
    abstol = 1e-7*np.ones(rho_init.shape)
    
    # rho = scpintg.odeint(MasterEquation, rho0, tspan, args=(H, ScaledCollapseOperators,  w0, num), rtol=reltol, atol=abstol)
                                                          
    solver = scpintg.ode(MasterEquation).set_integrator('vode', method='bdf', atol=abstol, rtol=reltol)
    solver.set_initial_value(rho_init, tspan[0])
    solver.set_f_params(H, ScaledCollapseOperators, w0, num)
    

    num_steps = np.size(tspan)
    delta_t = tspan[1]-tspan[0]

    rho = np.zeros((num_steps, np.size(rho0)), dtype=np.complex128)
    rho[0, :] = rho0
    
    t = np.zeros((num_steps,))
    t[0] = tspan[0]
    
    k = 1
    while k < num_steps: # solver.successful() and
        solver.integrate(solver.t + delta_t)
        
        # Store the results
        t[k] = solver.t
        rho[k, :] = solver.y[0:num**2] + 1j*solver.y[num**2:]
        k += 1

    return rho, t



def LindbladFun(x, y, rho):

    # Inputs are numpy arrays

    xdag = np.transpose(x)
    ydag = np.transpose(y)

    out = 1./2.*(np.dot(np.dot(x, rho), ydag) + np.dot(np.dot(y, rho), xdag) - \
      np.dot(np.dot(ydag, x), rho) - np.dot(np.dot(rho, xdag), y))

    return out



# Plank function
def PlankFun(w, T): 
    
    # Planck's Constants (J.s)
    h = 6.626068e-34
    
    # Reduced Planck's Constant
    hbar = h/2.0/pi
    
    # Boltzman constant
    kb = 1.3806503e-23 # (Joules/Kelvin)   
    
    if T==0:
        N = 0.0
    else:
        N = 1.0/(np.exp(hbar*w/kb/T) - 1.0)
    return N







