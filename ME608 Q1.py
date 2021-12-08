"""
This code is used to solve the first question of Fundamentals of Combustion (ME608) group project.
Developed by:
Anmoldeep Singh 180030002
"""
# importing relevant libraries
import matplotlib.pyplot as plt
import numpy as np

# Defining required constants
m = 2  # number of Carbon atoms in fuel
n = 6  # number of Hydrogen atoms in fuel
phi = 0.6  # equivalence ratio
Diameter = 0.03  # in m
T_i = 1000  # in K
P_i = 101325  # in Pa
m_dot = 0.2  # in kg/s
epsilon = 0.8  # quantifies the extent of reaction
del_Hc = 4 * (10 ** 7)  # in J/kg
Cp = 0.122 * (10 ** 4)  # in J/kg-K
Af = 4.338 * (10 ** 8)  # pre-exponential factor
x = 0.1
y = 1.65
E = 15098  # in K
a_stoic = m + (n / 4)
Area = np.pi * (Diameter ** 2) / 4  # in m2
mol_mass_O2 = 32.0  # in kg/kmol
mol_mass_N2 = 28.0  # in kg/kmol
mol_mass_Oxidiser = (mol_mass_O2 + 3.76 * mol_mass_N2) / 4.76  # in kg/kmol
mol_mass_Fuel = m * 12.0 + n * 1.0  # in kg/kmol
mol_mass_CO2 = 44.0  # in kg/kmol
mol_mass_H2O = 18.0  # in kg/kmol
mol_mass_Product = (m * mol_mass_CO2 + (n / 2) * mol_mass_H2O) / (m + (n / 2))  # in kg/kmol
nu = (mol_mass_O2 * a_stoic) / mol_mass_Fuel
Ru = 8.315  # in kJ/kmol-K

# Defining initial state
# at x = 0
m_mix_i = mol_mass_Fuel + 4.76 * (a_stoic * mol_mass_Oxidiser) / phi  # initial mass of the mixture
Y_Fuel_i = mol_mass_Fuel / m_mix_i  # initial mass fraction of fuel
Y_O2_i = (a_stoic * mol_mass_O2) / (phi * m_mix_i)  # initial mass fraction of oxygen
Y_Prod_i = 0  # initial mass fraction of products
Y_N2_i = 1 - Y_Fuel_i - Y_O2_i - Y_Prod_i  # initial mass fraction of nitrogen
mol_mass_mix = m_mix_i / (1 + 4.76 * a_stoic / phi)  # initial molecular mass of the mixture

# Defining final state
# at x = L
M_mix_f = ((1 - epsilon) * mol_mass_Fuel + epsilon * (m + n / 2) * mol_mass_Product + (  # final mass of the mixture
            (a_stoic / phi) - epsilon * (m + n / 2)) * mol_mass_O2 + ((3.76 * a_stoic * mol_mass_N2) / phi))
Y_Fuel_f = (1 - epsilon) * (mol_mass_Fuel / M_mix_f)  # final mass fraction of fuel
Y_O2_f = ((a_stoic / phi) - epsilon * (m + n / 2)) * mol_mass_O2 / M_mix_f  # final mass fraction of oxygen
Y_Prod_f = (epsilon * (m + n / 2) * mol_mass_Product) / M_mix_f  # final mass fraction of products
Y_N2_f = Y_N2_i  # final mass fraction of nitrogen


# function to calculate the volumetric mass flow rate of fuel
def m_dot_triple_dash(t, m_mix):
    mf = Af * (((P_i * m_mix) / (1000*Ru * t)) ** (x + y)) * np.exp(-E / t) * (
                (Y_Fuel_i - (Cp * (t - T_i)) / del_Hc) ** x) * ((Y_O2_i - (nu * Cp * (t - T_i)) / del_Hc) ** y)
    return mf


# function to calculate the mass of the mixture
def mass_mix(y_f, y_o2, y_prod, y_n2):
    M = (y_f / mol_mass_Fuel) + (y_o2 / mol_mass_O2) + (y_prod / mol_mass_Product) + (y_n2 / mol_mass_N2)
    M_mix = 1 / M
    return M_mix


# initialization
step = 0  # initial step
del_step = 0.0001  # step size
T_old = T_i
M_mix_old = mol_mass_mix
Y_Fuel_old = Y_Fuel_i
Y_O2_old = Y_O2_i
Y_N2_old = Y_N2_i
Y_Prod_old = Y_Prod_i

# declaring lists to store values
step_list = [0]
T_list = [T_old]
Mix_mass_list = [M_mix_old]
Y_Fuel_list = [Y_Fuel_old]
Y_O2_list = [Y_O2_old]
Y_N2_list = [Y_N2_old]
Y_Prod_list = [Y_Prod_old]

i = 0  # to track the number of iterations
while Y_Prod_old < Y_Prod_f:  # iterating until calculated mass fraction of product is greater than theoretical
    i += 1
    # calculate the new temperature
    T_new = T_old + (del_step * Area * del_Hc * m_dot_triple_dash(T_old, M_mix_old)) / (m_dot * Cp)

    # calculate the new fuel mass fraction
    Y_Fuel_new = Y_Fuel_i + (Cp * (T_i - T_new)) / del_Hc

    # calculate the new oxygen mass fraction
    Y_O2_new = Y_O2_i + (nu * Cp * (T_i - T_new)) / del_Hc

    # calculate the new product mass fraction
    Y_Prod_new = Y_Prod_i - ((1 + nu) * Cp * (T_i - T_new)) / del_Hc

    # calculate the new nitrogen mass fraction
    Y_N2_new = Y_N2_old
    # Y_N2_new = 1 - Y_Fuel_new - Y_O2_new - Y_Prod_new

    # calculate new mass of mixture
    M_mix_new = mass_mix(Y_Fuel_new, Y_O2_new, Y_Prod_new, Y_N2_new)

    # appending values in lists
    T_list.append(T_new)
    step_list.append(i * del_step)
    Y_Fuel_list.append(Y_Fuel_new)
    Y_O2_list.append(Y_O2_new)
    Y_Prod_list.append(Y_Prod_new)
    Y_N2_list.append(Y_N2_new)
    Mix_mass_list.append(M_mix_new)

    # Updating the values
    T_old = T_new
    Y_Fuel_old = Y_Fuel_new
    Y_O2_old = Y_O2_new
    Y_Prod_old = Y_Prod_new
    Y_N2_old = Y_N2_new
    M_mix_old = M_mix_new

print("Number of iterations: ", i, end='\n')  # number of iterations

# interpolating to calculate the required length
length_req = step_list[len(step_list)-2] + (step_list[len(step_list)-1] - step_list[len(step_list)-2])*((Y_Prod_f - Y_Prod_list[len(Y_Prod_list)-2])/(Y_Prod_list[len(Y_Prod_list)-1] - Y_Prod_list[len(Y_Prod_list)-2]))
print("The required length is: ", length_req)


# Plotting the data
try:
    # variation of product mass fraction with length of chamber
    plt.plot(step_list, Y_Prod_list, 'r', lw=2.5, label="Y_Pr")
    # plt.plot(step_list[len(step_list) - 1], Y_Prod_f, 'o')
    plt.xlabel('Length of chamber (m)', fontsize=20)
    plt.xticks(fontsize=15)
    plt.ylabel('Product mass fraction', fontsize=20)
    plt.yticks(fontsize=15)
    plt.title('Variation of Product Mass Fraction with Length', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid()
    plt.show()

    # variation of oxidiser mass fraction with length of chamber
    plt.plot(step_list, Y_O2_list, 'c', lw=2.5, label="Y_O2")
    # plt.plot(step_list[len(step_list) - 1], Y_O2_f, 'o')
    plt.xlabel('Length of chamber (m)', fontsize=20)
    plt.xticks(fontsize=15)
    plt.ylabel('Oxygen mass fraction', fontsize=20)
    plt.yticks(fontsize=15)
    plt.title('Variation of Oxygen Mass Fraction with Length', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid()
    plt.show()

    # variation of fuel mass fraction with length of chamber
    plt.plot(step_list, Y_Fuel_list, 'm', lw=2.5, label="Y_Fuel")
    # plt.plot(step_list[len(step_list) - 1], Y_Fuel_f, 'o')
    plt.xlabel('Length of chamber (m)', fontsize=20)
    plt.xticks(fontsize=15)
    plt.ylabel('Fuel mass fraction', fontsize=20)
    plt.yticks(fontsize=15)
    plt.title('Variation of Fuel Mass Fraction with Length', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid()
    plt.show()

    # variation of Temperature with length of chamber
    plt.plot(step_list, T_list, 'g', lw=2.5, label="T")
    plt.xlabel('Length of chamber (m)', fontsize=20)
    plt.xticks(fontsize=15)
    plt.ylabel('Temperature (K)', fontsize=20)
    plt.yticks(fontsize=15)
    plt.title('Variation of Temperature with Length', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid()
    plt.show()

    # variation of mixture mass with length of chamber
    plt.plot(step_list, Mix_mass_list, 'b', lw=2.5, label="Molecular mass of mixture")
    plt.xlabel('Length of chamber (m)', fontsize=20)
    plt.xticks(fontsize=15)
    plt.ylabel('Molecular mass of mixture (kg/kmol)', fontsize=20)
    plt.yticks(fontsize=15)
    plt.title('Variation of Molecular Mass of Mixture with Length', fontsize=20)
    plt.legend(fontsize=20)
    plt.grid()
    plt.show()
finally:
    pass

# to plot variation of all quantities with length of chamber
plt.plot(step_list, Y_Prod_list, 'r', lw=2.5, label="Y_Pr")
plt.plot(step_list, Y_O2_list, 'c', lw=2.5, label="Y_O2")
plt.plot(step_list, Y_Fuel_list, 'm', lw=2.5, label="Y_Fuel")
plt.plot(step_list, Y_N2_list, lw=2.5, label="Y_N2")
plt.plot(step_list, [T/1000 for T in T_list], 'g', lw=2.5, label="T/1000")  # scaling temperature down to fit on the graph
plt.plot(step_list, [M/20 for M in Mix_mass_list], 'b', lw=2.5, label="Molecular mass of mixture/20")  # scaling molecular mass of mixture down to fit on the graph
plt.xlabel('Length of chamber (m)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Quantities', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Quantities', fontsize=20)
plt.legend(fontsize=12)
plt.ylim([-10e-2, 2.5])
plt.margins(0)
plt.grid()
plt.show()
