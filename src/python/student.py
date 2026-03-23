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

    # If no output grid is given, make a basic one.
    if t_eval is None:
        t_eval = np.linspace(t0, tf, 200)

    t_eval = np.array(t_eval, dtype=float)

    # Number of state variables
    n = len(y0)

    # Number of time points
    m = len(t_eval)

    # Storage for solution
    Y = np.zeros((n, m), dtype=float)

    # Start at the initial value
    y = y0.copy()
    Y[:, 0] = y

    # Go from one requested time to the next using one RK4 step each time.
    for j in range(m - 1):
        t = t_eval[j]
        h = t_eval[j + 1] - t_eval[j]

        k1 = np.array(fun(t, y), dtype=float)
        k2 = np.array(fun(t + 0.5 * h, y + 0.5 * h * k1), dtype=float)
        k3 = np.array(fun(t + 0.5 * h, y + 0.5 * h * k2), dtype=float)
        k4 = np.array(fun(t + h, y + h * k3), dtype=float)

        y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        Y[:, j + 1] = y

    return SimpleSolution(t_eval, Y)


def solve_continuous_are(A, B, Q, R):
    """
    Solve the continuous-time algebraic Riccati equation

        A^T P + P A - P B R^{-1} B^T P + Q = 0

    using the Hamiltonian matrix method.

    This is a standard direct method from control theory.
    """

    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    Q = np.array(Q, dtype=float)
    R = np.array(R, dtype=float)

    n = A.shape[0]

    # Compute R^{-1}
    Rinv = np.linalg.inv(R)

    # Build Hamiltonian matrix
    H = np.block([
        [A, -B @ Rinv @ B.T],
        [-Q, -A.T]
    ])

    # Compute eigenvalues and eigenvectors
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

    # Compute P = U2 * U1^{-1}
    P = U2 @ np.linalg.inv(U1)

    # Remove tiny imaginary parts from roundoff
    P = np.real(P)

    # Make P symmetric, since the Riccati solution should be symmetric
    P = 0.5 * (P + P.T)

    return P
