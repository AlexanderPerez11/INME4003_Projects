import numpy as np

La=10
Lb = 15
Da = 0.1
D1 = 0.08
D2 = 0.06
eb = 0.00006
ea = 0.000045
Qt = 0.003
eta_m = 0.7
Za = 25
Zb = 5

rho = 998
mu = 0.001002
nu = mu/rho

Cv = 340
C_elbow = 30
g = 9.81

def f_haaland(Re,D,e):
    return 0.3086/((np.log10(6.9/Re+((e/D)/3.9)**(1.11)))**2)

def Reynolds(Q,D,nu):
    return 4*Q/(np.pi*D*nu)

def losses(f, Q, D, L):
    return 8*Q**2/(np.pi**2*g*D**4)*(f*L/D)

Q1 = np.linspace(0,Qt, 300)
Q2 = np.linspace(0,Qt, 300)

Re1 = np.zeros_like(Q1)
Re2 = np.zeros_like(Q1)
f1 = np.zeros_like(Q1)
f2 = np.zeros_like(Q1)
hf1 = np.zeros_like(Q1)
hf2 = np.zeros_like(Q1)
match_index = [0 for i in range(len(Q1))]
print(Q1)
for i in range(len(Q1)-1):
    Re1[i] = Reynolds(Q1[i+1],D1,nu)
    Re2[i] = Reynolds(Q2[i+1], D2, nu)
    f1[i] = f_haaland(Re1[i+1],D1,eb)
    f2[i] = f_haaland(Re2[i+1], D2, ea)
    hf1[i] = losses(f1[i+1], Q1[i+1],D1,Lb)
    hf2[i] = losses(f2[i+1], Q2[i+1], D2, Lb)

print(Re1)


for i in range(len(Q1)):

    plus_percent = hf1[i]+0.05*hf1[i]
    minus_percent = hf1[i] - 0.05 * hf1[i]
    match_index[i] = np.where(hf2==2*hf1[i])


# print(match_index)