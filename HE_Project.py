import numpy as np
import matplotlib.pyplot as plt
import pyromat as pm

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




steam = pm.get('mp.H2O')

steam_state_in = steam.state(x=1, p=1.01325)
steam_state_out = steam.state(x=0, p=1.01325)

T_h_in = steam_state_in["T"][0]
T_h_out = steam_state_out["T"][0]
h_steam_in = steam_state_in["h"][0]
h_steam_out = steam_state_out["h"][0]
m_dot_s = 3900*(1/3600)  # mass flow rate of steam kg/s
Q = m_dot_s*(h_steam_in-h_steam_out)
print(Q)