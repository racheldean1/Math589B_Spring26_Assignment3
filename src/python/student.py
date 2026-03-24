# 
# In this module implement these two functions:
# 1. solve_continuous_are
# 2. solve_ivp
#
# Make sure that they are compatible with their usage
# in modal_lqr.py.
#
# The Gradescope Autograder will call your implementation
# through functions:
#
# 1. simulate_closed_loop
# 2. simulate_open_loop

import numpy as np

class SimpleSolution:    #For copying the part of SciPy's solve_ivp result that we want. (modallqr only uses sol.t and sol.y)
    def __init__(self, t, y):
        self.t = t
        self.y = y

def solve_ivp(fun, t_span, y0, t_eval=None, rtol=1e-8, atol=1e-10, **kwargs):    # ODE solver
    t0, tf = t_span
    y0 = np.array(y0, dtype=float)
    if t_eval is None:
        t_eval = np.linspace(t0, tf, 200)

    t_eval = np.array(t_eval, dtype=float)

    
    n = len(y0)        # number of state variables
    m = len(t_eval)       # number of time points

    # storage for solution
    Y = np.zeros((n, m), dtype=float)

    # Starting at the initial value
    y = y0.copy()
    Y[:, 0] = y

    def rk4_step(t, y, h):        #adding a helper function
        k1 = np.array(fun(t, y), dtype=float)
        k2 = np.array(fun(t + 0.5 * h, y + 0.5 * h * k1), dtype=float)
        k3 = np.array(fun(t + 0.5 * h, y + 0.5 * h * k2), dtype=float)
        k4 = np.array(fun(t + h, y + h * k3), dtype=float)
        return y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    

    # Go from one time to the next
    for j in range(m - 1):
        t_left = t_eval[j]
        t_right = t_eval[j + 1]
        h_total = t_right - t_left
        nsub = 10            #number of substeps 
        h = h_total / nsub

        t = t_left
        for _ in range(nsub):
            y = rk4_step(t, y, h)
            t = t + h

        Y[:, j + 1] = y
    return SimpleSolution(t_eval, Y)

def solve_continuous_are(A, B, Q, R): # Solving A^T P + P A - P B R^{-1} B^T P + Q = 0 using Hamiltonian matrix method
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    Q = np.array(Q, dtype=float)
    R = np.array(R, dtype=float)
    n = A.shape[0]

    # compute R^{-1}
    Rinv = np.linalg.inv(R)

    # Build Hamiltonian matrix
    H = np.block([
        [A, -B @ Rinv @ B.T],
        [-Q, -A.T]
    ])

    # compute eigenvalues and eigenvectors
    eigvals, eigvecs = np.linalg.eig(H)

    # Keep the eigenvectors for eigenvalues with negative real part.
    stable_indices = []
    for i in range(len(eigvals)):
        if np.real(eigvals[i]) < 0:
            stable_indices.append(i)

    if len(stable_indices) != n:
        raise ValueError("Could not find the correct stable subspace")

    U = eigvecs[:, stable_indices]

    # Split U into top half and bottom half
    U1 = U[:n, :]
    U2 = U[n:, :]
    P = U2 @ np.linalg.inv(U1)        # Compute P = U2 * U1^{-1}
    P = np.real(P)          # Remove imaginary parts 
    P = 0.5 * (P + P.T)             # Make P symmetric, since the Riccati solution should be symmetric

    return P
