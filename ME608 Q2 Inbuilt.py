"""
This code is an alternate solution used to solve the second question of Fundamentals of Combustion (ME608) group project.
Developed by:
Anmoldeep Singh 180030002
"""
# importing relevant libraries
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import numpy as np

# Defining required constants
m = 3  # number of Carbon atoms in fuel
n = 8  # number of Hydrogen atoms in fuel
phi = 1.0  # equivalence ratio
Diameter = 0.08  # in m
T_1 = 298.0  # in K
P_i = 101325  # in Pa
m_dot = 0.1  # in kg/s
M_mix = 29.0
Cp = 1230.0  # in J/kg-K
del_Hc = 4.65 * (10 ** 7)  # in J/kg
T_ad = 2394.0  # in K
Af = 4.836 * (10 ** 9)  # pre-exponential factor
x = 0.1
y = 1.65
Ru = 8.3145  # in kJ/kmol-K
E = 15098 * Ru  # in K
Vol = np.pi * (Diameter ** 3) / 6  # in m3
a_stoic = m + (n / 4)
mass_O2 = 32.0  # in kg/kmol
mass_N2 = 28.0  # in kg/kmol
mass_Fuel = m * 12.0 + n * 1.0  # in kg/kmol
nu = (mass_O2 * a_stoic) / mass_Fuel

# storing all constant values
alpha = ((del_Hc * Vol) / Cp) * Af * (((P_i * M_mix) / (1000 * Ru)) ** (x + y)) * (mass_Fuel ** (1 - x)) * (mass_O2 ** (-y))
beta = (-E / Ru)
gamma = Cp / del_Hc

# Setting up initial conditions
Y_Fuel_1 = mass_Fuel / (mass_Fuel + (a_stoic / phi) * (mass_O2 + 3.76 * mass_N2))  # initial mass fraction of fuel
Y_O2_1 = (mass_O2 * (a_stoic / phi)) / (mass_Fuel + (a_stoic / phi) * (mass_O2 + 3.76 * mass_N2))  # initial mass fraction of oxygen
Y_N2_1 = 1 - Y_Fuel_1 - Y_O2_1  # initial mass fraction of nitrogen


# to calculate value of the function
def func(T):
    func_val = T[0] - (alpha / T[1]) * np.exp(beta / T[0]) * (T[0] ** (-(x + y))) * ((Y_Fuel_1 - gamma * (T[0] - T_1)) ** x) * (
                (Y_O2_1 - nu * gamma * (T[0] - T_1)) ** y) - T_1
    return [func_val, 0]


# initial guesses
T_cold_guess = 500.0  # in K, initial guess for estimating cold temperature
T_unstable_guess = 1400.0  # in K, initial guess for estimating unstable temperature
T_hot_guess = 2400.0  # in K, initial guess for estimating hot temperature

# declaring lists to store values
T_list_1 = [T_cold_guess]
T_list_2 = [T_unstable_guess]
T_list_3 = [T_hot_guess]
Y_Fuel_list = []
Y_O2_list = []
Y_N2_list = []

# list to store values of mass flow rate
mass = [0.1]

# storing different values of mass flow rate
for l in range(1, 10000):
    mass.append(round((mass[len(mass)-1] + 0.0001), 4))

# iterating until either all the mass flow rates are used to blowout condition is met
for m in mass:

    # calculating the temperatures
    T_cold = fsolve(func, [T_list_1[len(T_list_1) - 1], m])
    T_unstable = fsolve(func, [T_list_2[len(T_list_2) - 1], m])
    T_hot = fsolve(func, [T_list_3[len(T_list_3) - 1], m])

    # calculating the mass fractions
    Y_Fuel_new = Y_Fuel_1 - (Cp * (T_unstable[0] - T_1)) / del_Hc
    Y_O2_new = Y_O2_1 - (nu*Cp*(T_unstable[0]-T_1))/del_Hc
    Y_N2_new = 1 - Y_Fuel_new - Y_O2_new

    # un-comment and run to print the values
    # print(m, ":", round(T_cold[0], 4), round(T_unstable[0], 4), round(T_hot[0], 4), sep="  ")

    # checking for blowout condition
    # error at the end should not be worried about because that is the stopping criteria.
    if abs(T_unstable[0] - T_hot[0]) < 10e-3:
        print("Blow off limit reached for mass flow rate = ", m, "kg/s")
        break

    # appending values in lists
    T_list_1.append(round(T_cold[0], 4))
    T_list_2.append(round(T_unstable[0], 4))
    T_list_3.append(round(T_hot[0], 4))
    Y_Fuel_list.append(Y_Fuel_new)
    Y_O2_list.append(Y_O2_new)
    Y_N2_list.append(Y_N2_new)

# removing initial guess to avoid abrupt jumps in graphs
T_list_1.remove(T_list_1[0])
T_list_2.remove(T_list_2[0])
T_list_3.remove(T_list_3[0])

mass = mass[:len(T_list_1)]  # making the length of mass flow rate list same as temperature list

# Plotting the data
plt.plot(mass, T_list_1, 'c', lw=2.5, label="Cold Temperature")  # variation of cold temperature with mass flow rate
plt.plot(mass, T_list_2, 'y', lw=2.5, label="Unstable Temperature")  # variation of unstable temperature with mass flow rate
plt.plot(mass, T_list_3, 'r', lw=2.5, label="Hot Temperature")  # variation of hot temperature with mass flow rate
plt.xlabel('Mass Flow Rate (kg/s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Temperature (K)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Temperatures with Mass flow rate', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# plotting variation of mass fractions with mass flow rate
plt.plot(mass, Y_O2_list, lw=2.5, label="Y_O2")  # variation of oxygen mass fraction with mass flow rate
plt.plot(mass, Y_Fuel_list, lw=2.5, label="Y_Fuel")  # variation of fuel mass fraction with mass flow rate
plt.plot(mass, Y_N2_list, lw=2.5, label="Y_N2")  # variation of nitrogen mass fraction with mass flow rate
plt.xlabel('Mass Flow Rate (kg/s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Mass Fractions', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Mass fractions with Mass flow rate', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# to print the final exit conditions
print("Final fuel mass fraction: ", Y_Fuel_list[len(Y_Fuel_list)-1])
print("Final oxygen mass fraction: ", Y_O2_list[len(Y_O2_list)-1])
print("Final nitrogen mass fraction: ", Y_N2_list[len(Y_N2_list)-1])
print("Sum of all mass fractions: ", Y_Fuel_list[len(Y_Fuel_list)-1]+Y_O2_list[len(Y_O2_list)-1]+Y_N2_list[len(Y_N2_list)-1])
print("Final cold temperature: ", round(T_list_1[len(T_list_1)-1], 4), "K")
print("Final unstable temperature: ", round(T_list_2[len(T_list_2)-1], 4), "K")
print("Final  hot temperature: ", round(T_list_3[len(T_list_3)-1], 4), "K")