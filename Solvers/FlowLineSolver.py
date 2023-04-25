import numpy as np
L1 = 3000
L2 = 3000
D1 = 1
D2 = 0.67
Za = 100
Zb = 80

e1 = 0.001
e2 = 0.0001
v = 0.00003

Ws = 50
dz = Zb-Za
g = 32.174

def f1(x,v,D,e,L):
    a = f2(x,v,D,e)
    return  -Ws*g+dz+ (L/(2*g*D))*x**2*a

def f2(x,v,D,e):
    return 0.3086*(np.log10(6.9*v/D*x**(-1)+((e/D)/3.9)**(1.11)))**2


def df1dx(x,v,D,e,L):
    a = f2(x,v,D,e)
    b = df2dx(x,v,D,e)
    return (L/(2*g*D))*x*(2*a+x*b)

def df2dx(x,v,D,e):
    return -(0.116411* v *np.log(0.220759*(e/D)**1.11 + (6.9 *v)/(D*x)))/(x *(0.031994*D*x*(e/D)**1.11 + v))

def newton_taphson(f1, df1dx,f2, df2dx, x0, v,D,e,L,tol = 1e-7, max_iter = 1000000):
    x = x0
    iter = 0
    error  = 1
    while error > tol:
        fx = f1(x,v,D,e,L)
        dfx = df1dx(x,v,D,e,L)
        x -= fx/ dfx

        error  = abs(fx)

        iter +=1
        print(iter)
        if iter >= max_iter:
            raise ValueError(f"Failed to converge after {max_iter} iterations.")
    print(f"Reached solution after {iter} iterations")
    print(f'Error: {error}')
    return x

# for i in range(100):
#     T_guess = i+273.15
#     T_sol_my = newton_taphson(f, dfdx, T_guess)
Q1 = newton_taphson(f1,df1dx, f2,df2dx,6, v,D1,e1,L1)
print (f"Solution: {round(Q1,3)}")