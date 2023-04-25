import numpy as np
from scipy.optimize import fsolve
#Declare Parameters
Tinf = 27+ 273.15#Ambient Temp in K
Esun = 724 #Avg Sun Irradiance in W/m^2
alpha = 0.8 #Panel absorbtivity
hconv = 7.429 #convection coefficient in W/m^2k
emiss = 0.8 #Panel emissivity
boltz = 5.67e-8 #Steffan Boltzman Constant in W/m^2-k^4
# x is the surface temperature Ts
T_guess = 70  +273.15 #Initial guess of Ts in K

def tempsol(x):
    return -x + Tinf + (alpha*Esun - emiss*boltz*((x)**4 - (Tinf)**4))/hconv

T_sol = fsolve(tempsol,T_guess)
print(T_sol)
def f(x):
    #return q/(U*As*F)*np.log((Thi-x)/(Tho-Tci)) - (Thi-x - (Tho-Tci))
    return  -x + Tinf + (alpha*Esun - emiss*boltz*((x)**4 - (Tinf)**4))/hconv

def dfdx(x):
    #return q / (U * As * F) * (1/(Thi - x)) + 1
    return -1 -4*emiss*boltz*(x)**3


def newton_taphson(f, dfdx, x0, tol = 1e-7, max_iter = 1000000):
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

# for i in range(100):
#     T_guess = i+273.15
#     T_sol_my = newton_taphson(f, dfdx, T_guess)
T_sol_my = newton_taphson(f,dfdx, T_guess)
print (f"Solution: {round(T_sol_my,3)}")