import itertools
import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
import pandas as pd
import sys

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


def drawProgressBar(percent, barLen=20):
    # percent float from 0 to 1.
    sys.stdout.write("\r")
    sys.stdout.write("[{:<{}}] {:.0f}%".format("=" * int(barLen * percent), barLen, percent * 100))
    sys.stdout.flush()


def h_i_calculator(V, D):
    Re = (rho_c * V * D) / mu_c  # Internal Flow Reynolds Number
    if Re < 2300:
        Nu = 3.66  # laminar Nusselt for constant temperature cylinder
        h = (Nu * k_c) / D  # Convection coefficient [W/m^2-K]
        return h
    else:
        n = 0.4
        Nu = 0.023 * Re ** 0.8 * Pr_c ** n  # Turbulent Nusselt number
        h = (Nu * k_c) / D
        return h


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

    D_bundle = d_o * (N_tubes / K1) ** (1 / n1)  # Bundle Diameter

    return D_bundle


def calaculate_properties(Tout, Tin):
    Tm = (Tout + Tin) / 2  # Mean Temperature

    Tf = (Tm + T_h) / 2  # Film Temperature

    Cooling_Water = Fluid(FluidsList.Water).with_state(Input.temperature(Tm), Input.quality(0))

    C_p_c = Cooling_Water.specific_heat  # Specific heat of cooling water [J/kg-K]
    rho_c = Cooling_Water.density  # Density cooling water [kg/m^3]
    k_c = Cooling_Water.conductivity  # Conductivity of cooling water [W/m-k]
    mu_c = Cooling_Water.dynamic_viscosity  # Dynamic viscosity of cooling water [kg/m-s]
    Pr_c = Cooling_Water.prandtl  # Prandlt number of cooling water

    Steam_mean = Fluid(FluidsList.Water).with_state(Input.temperature(Tf), Input.quality(100))
    C_p_h = Steam_mean.specific_heat  # Specific heat of cooling water [J/kg-K]
    rho_h = Steam_mean.density  # Density cooling water [kg/m^3]
    k_h = Steam_mean.conductivity  # Conductivity of cooling water [W/m-k]
    mu_h = Steam_mean.dynamic_viscosity  # Dynamic viscosity of cooling water [kg/m-s]
    Pr_h = Steam_mean.prandtl  # Prandlt number of cooling water

    return [C_p_c, C_p_h, rho_c, rho_h, k_c, k_h, mu_c, mu_h, Pr_c, Pr_h]


# Calculate Steam Properties
P_steam = 101325  # Steam Inlet pressure [Pa]

Steam_in = Fluid(FluidsList.Water).with_state(Input.pressure(P_steam), Input.quality(100))
Steam_out = Fluid(FluidsList.Water).with_state(Input.pressure(P_steam), Input.quality(0))

T_h = Steam_in.temperature  # Constant Phase Change Temperature of steam [C]

h_fg_steam = Steam_in.enthalpy - Steam_out.enthalpy  # Vaporization Enthalpy of water [J/kg-K]
T_c_in = 40  # temperature of cooling water at the inlet [C]

m_dot_h = 3900 / 3600  # mass flow rate of steam [kg/s]
Q_req = m_dot_h * (h_fg_steam)  # Heat needed to condense the steam [W]

########################################################################################################################
U_o_guess = np.linspace(1100, 5600, 9)  # Assumed Overall HT Coefficient [W/m^2-k]

T_c_out_guess = np.linspace(T_c_in + 25, T_h - 10, 5)  # Outlet Temperature Array [C]
n_passes = np.array([2, 4, 6])  # Array of allowable tube passes
n_tubes = np.linspace(20, 200, 19)  # Array of tube numbers

material_k = np.array([237, 401, 25, 92, 60.5])  # Material conductivity [W/m-K]
material_name = ["Aluminum", "Copper", "Stainless Steel", "Nickel", "Carbon Steel"]
material_dict = {"Material": material_name, "Conductivity": material_k}
material_index = [i for i in range(len(material_name))]

tube_size = {
    "d_o_tubes": [1.05 * 0.0254, 1.315 * 0.0254, 1.66 * 0.0254, 1.9 * 0.0254, 2.375 * 0.0254, 2.875 * 0.0254,
                  3.5 * 0.0254,
                  4 * 0.0254, 3.5 * 0.0254],
    "d_i_tubes": [0.824 * 0.0254, 1.049 * 0.0254, 1.380 * 0.0254, 1.610 * 0.0254, 2.067 * 0.0254, 2.469 * 0.0254,
                  3.068 * 0.0254,
                  3.548 * 0.0254, 4.029 * 0.0254]}  # Allowable tube sizes [m]
tube_index = [i for i in range(len(tube_size["d_o_tubes"]))]
l_tube = np.linspace(4, 8, 4)  # Array of allowable tube lengths

# Combinations vector
tube_type_combs = list(
    itertools.product(n_passes, material_index, tube_index, l_tube, U_o_guess, T_c_out_guess, n_tubes))

# Arrays for storing data
combinations = len(tube_type_combs)
U = np.zeros((len(tube_type_combs)))
A = np.zeros((len(tube_type_combs)))
e = np.zeros((len(tube_type_combs)))
Pitch = np.zeros((len(tube_type_combs)))
N_tubes = np.zeros((len(tube_type_combs)))
V_tubes = np.zeros((len(tube_type_combs)))
V_shell = np.zeros((len(tube_type_combs)))
d_shell = np.zeros((len(tube_type_combs)))
conv_o = np.zeros((len(tube_type_combs)))
conv_i = np.zeros((len(tube_type_combs)))
d_outside = np.zeros((len(tube_type_combs)))
d_inside = np.zeros((len(tube_type_combs)))
N_p = np.zeros((len(tube_type_combs)))
L_tubes = np.zeros((len(tube_type_combs)))
k_names = np.zeros((len(tube_type_combs)))
m_dot_cold = np.zeros((len(tube_type_combs)))
dT_log_mean = np.zeros((len(tube_type_combs)))
A_req = np.zeros((len(tube_type_combs)))
U_o_g = np.zeros((len(tube_type_combs)))
T_c_o_g = np.zeros((len(tube_type_combs)))

count = 0

for array in tube_type_combs:
    # Extract Properties
    N_passes = array[0]

    d_o = tube_size["d_o_tubes"][array[2]]
    d_i = tube_size["d_i_tubes"][array[2]]
    P_t = 1.25 * d_o
    k_material = material_dict["Conductivity"][array[1]]

    L_t = array[3]
    U_guess = array[4]
    T_c_out = array[5]

    fluid_prop = calaculate_properties(T_c_out, T_c_in)
    C_p_c = fluid_prop[0]
    C_p_h = fluid_prop[1]
    rho_c = fluid_prop[2]
    rho_h = fluid_prop[3]
    k_c = fluid_prop[4]
    k_h = fluid_prop[5]
    mu_c = fluid_prop[6]
    mu_h = fluid_prop[7]
    Pr_c = fluid_prop[8]
    Pr_h = fluid_prop[9]

    m_dot_c = Q_req / (C_p_c * (T_c_out - T_c_in))
    dT_lm = (T_h - T_c_in - (T_h - T_c_out)) / np.log((T_h - T_c_in) / (T_h - T_c_out))
    A_required = Q_req / (U_guess * dT_lm)
    # Estimate Number of Tubes
    N_tubes[count] = array[6]

    # Tube Side Velocity
    V_tubes[count] = m_dot_c / ((N_tubes[count]) * (rho_c * 0.25 * np.pi * d_i ** 2))

    # Calculate Internal Convection coefficient
    h_i = h_i_calculator(V_tubes[count], d_i)

    # Bundle Diameter
    shell_clearance = 0.039
    D_bundle = bundle_diameter(d_o, N_passes, N_tubes[count])
    D_shell = (D_bundle + shell_clearance)
    d_shell[count] = D_shell
    I_baffle = 1.6 * D_shell
    V_shell[count] = m_dot_h / (rho_h * D_shell * I_baffle)
    N_r = int(2 * N_tubes[count] / 3)

    rho_l = Steam_out.density
    rho_v = Steam_in.density
    k_l = Steam_out.conductivity
    mu_l = Steam_out.dynamic_viscosity
    g = 9.81
    gamma_c = m_dot_h / (L_t * N_tubes[count])

    h_o = 0.95 * k_l * (rho_l * (rho_l - rho_v) * 9.81 / (mu_l * gamma_c)) ** (1 / 3) * (N_r) ** (-1 / 6)

    conv_i[count] = h_i
    conv_o[count] = h_o
    RoAo = 1 / (h_o)
    RiAo = (1 / (h_i)) * (d_o / d_i)
    RwAo = d_o * np.log(d_o / d_i) / (2 * k_material)

    U[count] = 1 / (RiAo + RoAo + RwAo)
    A[count] = d_o * np.pi * L_t * N_tubes[count]
    NTU = (U[count] * A[count]) / (C_p_c * m_dot_c)
    e[count] = 1 - np.exp(-NTU)
    Pitch[count] = P_t
    d_outside[count] = d_o
    d_inside[count] = d_i
    N_p[count] = N_passes
    L_tubes[count] = L_t
    k_names[count] = array[1]
    m_dot_cold[count] = m_dot_c
    dT_log_mean[count] = dT_lm
    A_req[count] = A_required
    U_o_g[count] = U_guess
    T_c_o_g[count] = T_c_out
    drawProgressBar((count / combinations))
    count += 1

data_U = pd.DataFrame(data=U, columns=["U"])

data_d_shell = pd.DataFrame(data=d_shell, columns=["d_shell"])
data_d_o = pd.DataFrame(data=d_outside, columns=["d_o"])
data_d_i = pd.DataFrame(data=d_inside, columns=["d_i"])
data_V_tubes = pd.DataFrame(data=V_tubes, columns=["V_tubes"])
data_V_shell = pd.DataFrame(data=V_shell, columns=["V_shell"])
data_conv_o = pd.DataFrame(data=conv_o, columns=["h_o"])
data_conv_i = pd.DataFrame(data=conv_i, columns=["h_i"])
data_N_tubes = pd.DataFrame(data=N_tubes, columns=["N_t"])
data_area = pd.DataFrame(data=A, columns=["A"])
data_pitch = pd.DataFrame(data=Pitch, columns=["P_t"])
data_effectiveness = pd.DataFrame(data=e, columns=["e"])
data_L_t = pd.DataFrame(data=L_tubes, columns=["L_t"])
data_N_p = pd.DataFrame(data=N_p, columns=["N_p"])
data_names = pd.DataFrame(data=k_names, columns=["Material"])

data_A_req = pd.DataFrame(data=A_req, columns=["A_req"])
data_U_o_guess = pd.DataFrame(data=U_o_g, columns=["U_o_guess"])
data_T_o_guess = pd.DataFrame(data=T_c_o_g, columns=["T_o"])
data_dT_lm = pd.DataFrame(data=dT_log_mean, columns=["dT_lm"])
data_m_dot_c = pd.DataFrame(data=m_dot_cold, columns=["m_dot_c"])
data_HE = pd.concat(
    [data_names, data_N_p, data_L_t, data_pitch, data_N_tubes, data_d_o, data_d_i, data_area, data_d_shell,
     data_V_tubes, data_V_shell,
     data_conv_o, data_conv_i,
     data_U, data_effectiveness, data_A_req, data_U_o_guess, data_T_o_guess, data_dT_lm, data_m_dot_c],
    axis=1)
########################################################################################################################
filtered_index = []
for i in range(len(data_HE["d_shell"])):
    if 1100 < data_HE["U"][i] < 5600 and 0.9 < data_HE["V_tubes"][i] < 2.5 and data_HE["e"][i] >= 0.6 and 10 < \
            data_HE["V_shell"][i] < 30 and 1 / 15 < data_HE["D_s"][i] / data_HE["L_t"][i] < 1 / 5:  #
        filtered_index.append(i)

filtered_data = data_HE.loc[filtered_index]
# filtered_data.to_csv("filtered_HE_Data.csv")
########################################################################################################################
plt.figure(1)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["V_shell"], marker='.', label="shell")
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["V_tubes"], marker='.', label="tubes")
plt.title("Tube and Shell Velocities vs. Combinations")
plt.ylabel("Velocity (m/s)")
plt.xlabel("Combination")
plt.legend()

plt.figure(2)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["h_o"], marker='.', label="h_o")
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["h_i"], marker='.', label="h_i")
plt.title("Tube and Shell Convection Coefficients vs. Combinations")
plt.ylabel(r"$(\frac{W}{m^2-K})$")
plt.xlabel("Combination")
plt.legend()

plt.figure(3)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["U"], marker='.', label="U")
plt.title("Overall HT Coefficient vs. Combinations")
plt.ylabel(r"$(\frac{W}{m^2-K})$")
plt.xlabel("Combination")
plt.legend()

plt.figure(4)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["d_shell"], marker='.', label="d_shell")
plt.title("D_s vs. Combinations")
plt.ylabel(r"$(m)$")
plt.xlabel("Combination")
plt.legend()

plt.figure(5)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["N_t"], marker='.', label="N_t")
plt.title("Tube number vs. Combinations")
plt.ylabel("Tube Number")
plt.xlabel("Combination")
plt.legend()

plt.figure(6)
plt.scatter([i for i in range(len(filtered_data["U"]))], filtered_data["L_t"], marker='.', label="L_t")
plt.title("Tube Length vs. Combinations")
plt.ylabel("m")
plt.xlabel("Combination")
plt.legend()

plt.show()

print(f"Phase Change Water Temperature: T_h = {round(T_h, 2)} [C]")
print(f"Vaporization Enthalpy of Water: h_fg = {round(h_fg_steam, 2)} [J/kg-K]")
print(f"Inlet Water Temperature: T_h = {round(T_c_in, 2)} [C]")
print(f"Selected Outlet Water Temperature: T_h = {round(T_c_out, 2)} [C]")
print(f"Steam Mass Flow Rate: m_dot = {round(m_dot_h, 2)} [kg/s]")
print(f'Required Heat Rate: Q_req = {round(Q_req, 2)} [W]')
print(f"Cooling Water Mass Flow Rate: m_dot = {round(m_dot_c, 2)} [kg/s]")
print(f'Cooling Water Specific Heat: C_p = {round(C_p_c, 2)} [J/kg-K]')
print(f'Cooling Water Density: rho = {round(rho_c, 2)} [kg/m^3]')
print(f'Cooling Water Conductivity: k = {round(k_c, 2)} [W/m-k]')
print(f'Cooling Water viscosity: mu = {round(mu_c, 5)} [kg/m-s]')
print(f'Cooling Water Prandlt-#: Pr = {round(Pr_c, 2)}')

print(f'Steam Specific Heat: C_p = {round(C_p_h, 2)} [J/kg-K]')
print(f'Steam Density: rho = {round(rho_h, 2)} [kg/m^3]')
print(f'Steam Conductivity: k = {round(k_h, 2)} [W/m-k]')
print(f'Steam viscosity: mu = {round(mu_h, 5)} [kg/m-s]')
print(f'Steam Prandlt-#: Pr = {round(Pr_h, 2)}')

print(f'Log Mean Temperature = {round(dT_lm, 2)} [C]')
print(f'Required Heat Rate: Q_req = {round(Q_req, 2)} [W]')
print(f'Required Area: A_req = {round(A_required, 2)} [m^2]')
