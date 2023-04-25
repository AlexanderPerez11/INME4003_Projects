import numpy as np
from scipy.optimize import fsolve
q = 1064981.925
U = 41.27
As  = 225
F = 0.98
Thi = 250
Tho = 200
Tci = 85

Tco = 100
tol = 1e-3
error = 1

def tempsol(Tco, Thi, Tci, Tho):
    return q / (U * As * F) * np.log((Thi - Tco) / (Tho - Tci)) - ( Thi - Tco - (Tho-Tci))

T_guess = -5000

T_sol = 2
)
print(T_sol)

def f(x):
    return q/(U*As*F)*np.log((Thi-x)/(Tho-Tci)) - (Thi-x - (Tho-Tci))

def dfdx(x):
    return q / (U * As * F) * (1/(Thi - x)) + 1


def newton_taphson(f, dfdx, x0, tol = 1e-7, max_iter = 10000):
    x = x0
    iter = 0
    error  = 1
    while error > tol:
        fx = f(x)
        dfx = dfdx(x)
        x -= fx/ dfx

        error  = abs(fx)

        iter +=1
        if iter >= max_iter:
            raise ValueError(f"Failed to converge after {max_iter} iterations.")
    print(f"Reached solution after {iter} iterations")
    print(f'Error: {error}')
    return x





T_sol_my = newton_taphson(f,dfdx, T_guess)
print (f"Solution: {round(T_sol_my,3)}")