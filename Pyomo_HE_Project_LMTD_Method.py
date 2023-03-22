import itertools
import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
import pandas as pd
import warnings

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


def h_o_calculator(V, D, S_t, S_l, N,D_e, aligned):
    if aligned:
        V_max = V * S_t / (S_t - D)
        if N < 20:
            N_L = [0.70, 0.80, 0.86, 0.90, 0.92, 0.935, 0.95, 0.957, 0.963, 0.97, 0.973, 0.977, 0.98, 0.983, 0.987,
                   0.99,
                   0.9925, 0.995, 0.9975, 1]
            C_2 = N_L[int(N) - 1]

        else:
            C_2 = 1
            print("Hi")
        Re = (V_max * rho_h * D) / mu_h
        if 10 < Re <= 100:
            C = 0.8
            m = 0.4
        elif 100 < Re <= 1000:
            Nu = 0.3 + (0.62 * Re ** (1 / 2) * Pr_h ** (1 / 3)) / ((1 + (0.4 / Pr_h) ** (2 / 3)) ** (1 / 4)) * (
                    1 + (Re / 282000) ** (5 / 8)) ** (4 / 5)
            h = (Nu * k_h) / D
            return h
        elif 1000 < Re <= 2 * 10 ** 5:
            C = 0.27
            m = 0.63
        elif 2 * 10 ** 5 < Re <= 2 * 10 ** 6:
            C = 0.021
            m = 0.84
        else:
            raise Exception("Something Aint Right")
    else:
        S_d = np.sqrt(S_l ** 2 + (S_t / 2) ** 2)
        if 2 * (S_d - D) >= (S_t - D):
            V_max = (S_t * V) / (S_t - D)
        else:
            V_max = (S_t * V) / (2 * (S_d - D))

        Re = (V_max * rho_h * D_e) / mu_h
        if N < 20:
            N_L = [0.64, 0.76, 0.84, 0.89, 0.92, 0.935, 0.95, 0.957, 0.963, 0.97, 0.973, 0.977, 0.98, 0.983, 0.987,
                   0.99,
                   0.9925, 0.995, 0.9975, 1]
            C_2 = N_L[int(N) - 1]
        else:
            C_2 = 1
            # print("hi")

        if 10 < Re <= 100:
            C = 0.90
            m = 0.4
        elif 100 < Re <= 1000:
            Nu = 0.3 + (0.62 * Re ** (1 / 2) * Pr_h ** (1 / 3)) / ((1 + (0.4 / Pr_h) ** (2 / 3)) ** (1 / 4)) * (
                    1 + (Re / 282000) ** (5 / 8)) ** (4 / 5)
            h = (Nu * k_h) / D
            return h
        elif 1000 < Re <= 2 * 10 ** 5:
            C = 0.35 * (S_t / S_l) ** 0.2
            m = 0.60
        elif 2 * 10 ** 5 < Re <= 2 * 10 ** 6:
            C = 0.022
            m = 0.84
        else:
            warnings.warn(f"Reynolds Number exceeds 2x10^6, Re= {Re}. Setting h to 0")
            return 1

    Nu = C_2 * (C * (Re ** m) * (Pr_h ** 0.36) * (Pr_h / Pr_s) ** (1 / 4))
    h = (Nu * k_h) / D
    return h


def h_i_calculator(V, D):
    # V = m_dot / (N* rho_h * np.pi * 0.25 * (D) ** 2)
    Re = (rho_c * V * D) / mu_c
    if Re < 2300:
        Nu = 3.66
        h = (Nu * k_c) / D
        return h
    else:
        if 90 > T_m:
            n = 0.4
        else:
            n = 0.3
        Nu = 0.023 * Re ** 0.8 * Pr_c ** n
        h = (Nu * k_c) / D
        return h


def steam_properties():
    pass


def bundle_diameter(d_o, N_passes, N_tubes):
    K = [0.319, 0.249, 0.175, 0.0743, 0.0365]
    n = [2.142, 2.207, 2.285, 2.499, 2.675]
    if N_passes == 1:
        K1 = K[0]
        n1 = n[0]

    if N_passes == 2:
        K1 = K[1]
        n1 = n[1]

    if N_passes == 4:
        K1 = K[2]
        n1 = n[2]

    if N_passes == 6:
        K1 = K[3]
        n1 = n[3]

    if N_passes == 8:
        K1 = K[4]
        n1 = n[4]

    D_bundle = d_o * (N_tubes / K1) ** (1 / n1)

    return D_bundle


# Calculate Steam Properties
P_steam = 101325  # Steam Inlet pressure [Pa]

Steam_in = Fluid(FluidsList.Water).with_state(Input.pressure(P_steam), Input.quality(0))
Steam_out = Fluid(FluidsList.Water).with_state(Input.pressure(P_steam), Input.quality(100))

T_h = Steam_in.temperature  # Constnat Phase Change Temperature of steam [C]

h_fg_steam = Steam_out.enthalpy - Steam_in.enthalpy  # Vaporization Enthalpy of water [J/kg-K]
print(Steam_out.enthalpy, Steam_in.enthalpy)
T_c_in = 40  # temperature of cooling water at the inlet [C]
T_c_out = 60  # temperature of cooling water at the outlet [C]

m_dot_h = 3900 / 3600 # mass flow rate of steam [kg/s]
print(h_fg_steam)
Q_req = m_dot_h * (h_fg_steam)  # Heat needed to condense the steam [W]

########################################################################################################################
# Calculate Cooling Water Properties at T_m
T_m = (T_c_in + T_c_out) / 2  # Mean Cooling Water Temperature

Cooling_Water = Fluid(FluidsList.Water).with_state(Input.temperature(T_m), Input.quality(0))

C_p_c = Cooling_Water.specific_heat  # Specific heat of cooling water [J/kg-K]
rho_c = Cooling_Water.density  # Density cooling water [kg/m^3]
k_c = Cooling_Water.conductivity  # Conductivity of cooling water [W/m-k]
mu_c = Cooling_Water.dynamic_viscosity  # Dynamic viscosity of cooling water [kg/m-s]
Pr_c = Cooling_Water.prandtl  # Prandlt number of cooling water

# print(f'Cooling Water Specific Heat: C_p = {round(C_p_c, 2)} [J/kg-K]')
# print(f'Cooling Water Density: rho = {round(rho_c, 2)} [kg/m^3]')
# print(f'Cooling Water Conductivity: k = {round(k_c, 2)} [W/m-k]')
# print(f'Cooling Water viscosity: mu = {round(mu_c, 5)} [kg/m-s]')
# print(f'Cooling Water Prandlt-#: Pr = {round(Pr_c, 2)}')

########################################################################################################################
# Fluid Properties for Convection Coefficient Calculation
# Properties at T_f for external tube bank flow T_f = (T_h + T_s)/2

T_s = (T_c_in + T_c_out) / 2
T_f = (T_h + T_s) / 2

Steam_mean = Fluid(FluidsList.Water).with_state(Input.temperature(T_f), Input.quality(1))
C_p_h = Steam_mean.specific_heat  # Specific heat of cooling water [J/kg-K]
rho_h = Steam_mean.density  # Density cooling water [kg/m^3]
k_h = Steam_mean.conductivity  # Conductivity of cooling water [W/m-k]
mu_h = Steam_mean.dynamic_viscosity  # Dynamic viscosity of cooling water [kg/m-s]
Pr_h = Steam_mean.prandtl  # Prandlt number of cooling water

Steam_surface_temp = Fluid(FluidsList.Water).with_state(Input.temperature(T_s), Input.quality(1))
Pr_s = Steam_surface_temp.prandtl
########################################################################################################################
# Calculate water mass flow rate
m_dot_c = Q_req / (C_p_c * (T_c_out - T_c_in))  # mass flow rate of cooling water [kg/s]
k_steel = 401  # Conductivity of copper [W/m-k]
dT_lm = (T_h - T_c_in - (T_h - T_c_out)) / np.log((T_h - T_c_in) / (T_h - T_c_out))  # Log Mean Temperature [C]
print(dT_lm)
########################################################################################################################
"""
Dtermine Required Surface Area estimated with outside tube diameter
"""
U_o_guess = 1100 # Assumed Overall HT Coefficient [W/m^2-k]
A_required = Q_req / (U_o_guess * dT_lm)

# n_tubes_l = np.linspace(1, 20, 20)  # Array of allowable number of tube rows
# n_tubes_t = np.linspace(1, 20, 20)
# baffle_spacing = np.linspace(0.3, 0.6, 10)  # Baffle Spacing percentage]

# tube_size = {"d_o_tubes": [1.05 * 0.0254, 1.315 * 0.0254, 1.], "d_i_tubes": [0.884 * 0.0254, 1.097 * 0.0254],
#              "P_t": [0.03667, 0.0459]}  # Allowable tube sizes [m]


n_passes = np.array([1, 2, 4, 6, 8])  # Array of allowable tube passes
material_k = np.array([237, 401, 25, 92, 60.5])
material_index = ["Aluminum", "Copper", "Stainless Steel", "Nickel", "Carbon Steel"]
tube_size = {"d_o_tubes": [1 * 0.0254, 1.25 * 0.0254, 1.5 * 0.0254, 1.75 * 0.0254],
             "d_i_tubes": [0.8 * 0.0254, 1.05 * 0.0254, 1.3 * 0.0254, 1.55 * 0.0254]}  # Allowable tube sizes [m]
tube_index = [i for i in range(len(tube_size["d_o_tubes"]))]
# l_tube = np.array([2.4384, 3.6576, 4.572, 6.096])  # Length of tubes [m]
l_tube = np.linspace(2, 6, 10)
tube_type_combs = list(itertools.product(n_passes, material_k, tube_index, l_tube))

U = np.zeros((len(tube_type_combs)))
N_tubes = np.zeros((len(tube_type_combs)))
# Q = np.zeros((len(tube_type_combs)))
V_tubes = np.zeros((len(tube_type_combs)))
V_shell = np.zeros((len(tube_type_combs)))
d_shell = np.zeros((len(tube_type_combs)))
conv_o = np.zeros((len(tube_type_combs)))
conv_i = np.zeros((len(tube_type_combs)))
count = 0

for array in tube_type_combs:
    N_passes = array[0]

    d_o = tube_size["d_o_tubes"][array[2]]
    d_i = tube_size["d_i_tubes"][array[2]]
    P_t = 1.25 * d_o
    k_material = array[1]

    L_t = array[3]
    # Estimate Number of Tubes
    N_tubes[count] = int(A_required / (d_o * L_t *np.pi))

    # Tube Side Velocity
    V_tubes[count] = m_dot_c / ((N_tubes[count] / N_passes) * (rho_c * 0.25 * np.pi * d_i ** 2))

    # Calculate Internal Convection coefficient
    h_i = h_i_calculator(V_tubes[count], d_i)

    # Bundle Diameter
    shell_clearance = 0
    D_bundle = bundle_diameter(d_o, N_passes, N_tubes[count])
    D_shell = D_bundle + shell_clearance
    d_shell[count] = D_shell
    I_baffle = D_shell
    # Cross Flow Area for tube banks
    A_s = (P_t - d_o) * D_shell * D_bundle / P_t
    G_shell = m_dot_h / A_s
    V_shell[count] = G_shell / rho_h
    d_e = 1.10 * (P_t ** 2 - 0.917 * d_o ** 2)
    N_r = int(2 * N_tubes[count] / 3)

    h_o = h_o_calculator(V_shell[count], d_o, P_t, P_t, N_r,d_e, False)

    conv_i[count] = h_i
    conv_o[count] = h_o
    RoAo = 1 / (h_o)
    RiAo = 1 / (h_i)
    RwAo = d_o * np.log(d_o / d_i) / (2 * k_material)

    U[count] = (RiAo + RoAo + RwAo) ** (-1)
    # baffle_percentage = array[5]
    # d_shell[count] = N_l * S_t * N_t

    # V_tubes[count] = m_dot_c / (N_t*N_l*rho_c * np.pi * 0.25 * d_i ** 2)
    # V_shell[count] = m_dot_h / (rho_c * d_shell[count] ** 2 * baffle_percentage)
    # h_o = h_o_calculator(V_shell[count], d_o, S_t, S_l, N_l, aligned=False)
    # h_i = h_i_calculator(V_tubes[count], d_i)
    # conv_i[count] = h_i
    # conv_o[count] = h_o
    # Ro = 1 / (N_t * N_l * d_o * L_t * h_o * np.pi)
    # Ri = 1 / (N_t * N_l * d_i * L_t * h_i * np.pi)
    # Rw = np.log(d_o / d_i) / (2 * np.pi * k_steel * L_t)
    # Ri = 1 *d_o/ (h_i*d_i)
    # Ro = 1 / (h_o)
    # Rw = np.log(d_o / d_i)*d_o / (2 * k_steel)
    # U[count] = (Ri + Ro + Rw) ** (-1)
    # A = N_t*N_l*d_o*np.pi*L_t
    # Q[count] = U[count] * dT_lm * A
    count += 1

data_U = pd.DataFrame(data=U, columns=["U"])
# data_Q = pd.DataFrame(data=Q, columns=["Q"])
data_DIM = pd.DataFrame(tube_type_combs, columns=["N_pass", "Material", "Tube_dim", "L_t"])
data_d_shell = pd.DataFrame(data=d_shell, columns=["d_shell"])
data_V_tubes = pd.DataFrame(data=V_tubes, columns=["V_tubes"])
data_V_shell = pd.DataFrame(data=V_shell, columns=["V_shell"])
data_conv_o = pd.DataFrame(data=conv_o, columns=["h_o"])
data_conv_i = pd.DataFrame(data=conv_i, columns=["h_i"])
data_N_tubes = pd.DataFrame(data=N_tubes, columns=["N_t"])

data_HE = pd.concat([data_U,data_DIM, data_N_tubes, data_d_shell, data_V_tubes, data_V_shell, data_conv_o, data_conv_i],
                    axis=1)
########################################################################################################################
filtered_index = []
for i in range(len(data_HE["d_shell"])):
    if 1100 < data_HE["U"][i] < 5600 and 0.9 <data_HE["V_tubes"][i]< 2.5 and 0.6 < data_HE["V_shell"][i] < 1.5:
        filtered_index.append(i)

filtered_data = data_HE.loc[filtered_index]
# ########################################################################################################################
plt.figure(1)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["V_shell"], marker='.', label="shell")
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["V_tubes"], marker='.', label="tubes")
plt.title("Tube and Shell Velocities vs. Combinations")
plt.ylabel("Velocity (m/s)")
plt.xlabel("Combination")
plt.legend()

plt.figure()
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["h_o"], marker='.', label="h_o")
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["h_i"], marker='.', label="h_i")
plt.title("Tube and Shell Convection Coefficients vs. Combinations")
plt.ylabel(r"$(\frac{W}{m^2-K})$")
plt.xlabel("Combination")
plt.legend()

plt.figure(3)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["d_shell"], marker='.', label="d_shell")
plt.title("Shell Diameter vs. Combinations")
plt.ylabel("r(m)")
plt.xlabel("Combination")
plt.legend()

plt.figure(4)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["N_t"], marker='.', label="d_shell")
plt.title("Tube number vs. Combinations")
plt.ylabel("r(m)")
plt.xlabel("Combination")
plt.legend()

plt.show()

print(f"Phase Change Water Temperature: T_h = {round(T_h, 2)} [C]")
print(f"Vaporization Enthalpy of Water: h_fg = {round(h_fg_steam, 2)} [J/kg-K]")
print(f"Inlet Water Temperature: T_h = {round(T_c_in, 2)} [C]")
print(f"Selected Outlet Water Temperature: T_h = {round(T_c_out, 2)} [C]")
print(f"Steam Mass Flow Rate: m_dot = {round(m_dot_h, 2)} [kg/s]")
print(f'Required Heat Rate: Q_req = {round(Q_req, 2)} [W]')
