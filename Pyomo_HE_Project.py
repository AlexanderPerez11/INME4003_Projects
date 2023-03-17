import itertools

import numpy as np
import matplotlib.pyplot as plt
import pyromat as pm
from pyomo.environ import *

"""
A shell and tube heat exchanger is to be designed to condense 3900 kg/hr of steam at 1 atm. The
steam may be assumed to enter the shell side as saturated vapor and leave as saturated liquid.
The overall heat transfer coefficient must be within the range of 1100-5600 𝑊
𝑚2𝐾
. Water is
employed as the cooling fluid in the tubes with an inlet temperature of 40℃.
Determine a suitable design that specifies:
(a) Outlet water temperature (℃)
(b) Number of shell passes
(c) Number of tube passes
(d) Tube diameter and length in each tube pass
(e) Mass flow rate of water (kg/hr)
(f) Shell inside diameter
(g) Condenser’s effectiveness
Your design must comply with standard parameters for this type of heat exchanger.
"""

steam = pm.get('mp.H2O') # Extract all data for water

steam_state_in = steam.state(x=1, p=1.01325) # get steam data at vapor saturation and atmospheric pressure
steam_state_out = steam.state(x=0, p=1.01325) # get steam data at liquid saturation and atmospheric pressure

T_h_in = steam_state_in["T"][0] # temperature of steam at the inlet [K]
T_h_out = steam_state_out["T"][0] # temperature of the steam at the outlet [K]
h_steam_in = steam_state_in["h"][0] # enthalpy of the steam at the inlet [kJ/kg]
h_steam_out = steam_state_out["h"][0]  # enthalpy of the steam at the outlet [kJ/kg]

T_c_in = 40 + 273.15 # temperature of cooling water at the inlet [K]

m_dot_s = 3900*(1/3600)  # mass flow rate of steam [kg/s]

Q_required = m_dot_s*(h_steam_in-h_steam_out) # Heat needed to condense the steam [W]
# print(f'Q_required = {round(Q_required,3)} W')
# print(steam_state_in)
# print(steam_state_out)

n_passes = np.linspace(1,10,10)
n_tubes = np.linspace(1,100,100)
d_o_tubes = np.array([15.88,19.05,25.4,32,39])
l_tube = np.array([2.4384,3.6576,4.572,6.096])
U_overall = np.linspace(1100,5600,45)

c_p_water = 4180 # water specific heat [J/kg-K]
rho_water  = 1000
n_sols = len(n_passes)* len(n_tubes)* len(d_o_tubes)* len(l_tube) * len(U_overall)

combinations = list(itertools.product(n_passes,n_tubes,d_o_tubes,l_tube,U_overall))
UA = np.zeros((len(combinations)))
NTU = np.zeros((len(combinations)))
effec = np.zeros((len(combinations)))
m_dot_w = np.zeros((len(combinations)))
v_water = np.zeros((len(combinations)))
T_c_out = np.zeros((len(combinations)))
print(UA.size)
print(len(combinations))
print(len(combinations[0]))
count = 0
for array in combinations:
    UA[count] = array[0]*array[1]*array[2]*array[3]*array[4]
    NTU[count] = UA[count]/c_p_water
    effec[count] = 1-np.exp(-NTU[count])
    m_dot_w[count] = Q_required/(effec[count]*c_p_water*(T_h_in - T_c_in))
    v_water[count] = m_dot_w[count]/(rho_water*array[0]*array[1]*array[2]*array[3])
    T_c_out = T_c_in+Q_required/(m_dot_w[count]*c_p_water)
    count+=1


variables_dict = {"UA":UA, "NTU":NTU, "effec":effec,"m_dot_water":m_dot_w,"T_c_out":T_c_out}
print(variables_dict["effec"])
print("Done!")

model = ConcreteModel()

model.Q = Var(bounds=())