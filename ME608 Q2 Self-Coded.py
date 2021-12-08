"""
This code is used to solve the second question of Fundamentals of Combustion (ME608) group project.
Developed by:
Anmoldeep Singh 180030002
"""
# importing relevant libraries
import matplotlib.pyplot as plt
import numpy as np

# Defining required constants
m = 3  # number of Carbon atoms in fuel
n = 8  # number of Hydrogen atoms in fuel
phi = 1.0  # equivalence ratio
Diameter = 0.08  # in m
T_1 = 298.0  # in K
P_1 = 101325  # in Pa
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
alpha = ((del_Hc * Vol) / Cp) * Af * (((P_1 * M_mix) / (1000 * Ru)) ** (x + y)) * (mass_Fuel ** (1 - x)) * (mass_O2 ** (-y))
beta = (-E / Ru)
gamma = Cp / del_Hc

# Setting up initial conditions
Y_Fuel_1 = mass_Fuel / (mass_Fuel + (a_stoic / phi) * (mass_O2 + 3.76 * mass_N2))  # initial mass fraction of fuel
Y_O2_1 = (mass_O2 * (a_stoic / phi)) / (mass_Fuel + (a_stoic / phi) * (mass_O2 + 3.76 * mass_N2))  # initial mass fraction of oxygen
Y_N2_1 = 1 - Y_Fuel_1 - Y_O2_1  # initial mass fraction of nitrogen


# to calculate value of the function
def func(T, m1):
    func_val = T - (alpha / m1) * np.exp(beta / T) * (T ** (-(x + y))) * ((Y_Fuel_1 - gamma * (T - T_1)) ** x) * (
                (Y_O2_1 - nu * gamma * (T - T_1)) ** y) - T_1
    return func_val


# to calculate value of the derivative of the function
def func_derivative(T, m2):
    func_derivative_val = 1 + (alpha/m2) * np.exp(beta / T) * (T ** (-(x + y))) * ((Y_Fuel_1 - gamma * (T - T_1)) ** x) * (
                (Y_O2_1 - nu * gamma * (T - T_1)) ** y) * (
                                 (beta / (T ** 2)) + ((x + y) / T) + ((x * gamma) / (Y_Fuel_1 - gamma * (T - T_1))) + (
                                     (y * nu * gamma) / (Y_O2_1 - nu * gamma * (T - T_1))))
    return func_derivative_val


# solving for temperature using Newton-Raphson method
def sol_for_T(t, m_):
    h = func(t, m_) / func_derivative(t, m_)
    while abs(h) >= 0.000001:
        h = func(t, m_) / func_derivative(t, m_)
        t = t - h
    return t


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
# change the iteration count (100000) and step size (0.0001) to get smaller observation points
for l in range(1, 100000):
    mass.append(round((mass[len(mass)-1] + 0.0001), 4))

# iterating until either all the mass flow rates are used to blowout condition is met
for m in mass:

    # calculating the temperatures
    T_cold = sol_for_T(T_list_1[len(T_list_1)-1], m)
    T_unstable = sol_for_T(T_list_2[len(T_list_2)-1], m)
    T_hot = sol_for_T(T_list_3[len(T_list_3)-1], m)

    # un-comment the line number 100 and run to print the values
    print(m, ":", round(T_cold, 4), round(T_unstable, 4), round(T_hot, 4), sep="  ")

    # calculating the mass fractions
    Y_Fuel_2 = Y_Fuel_1 - (Cp*(T_unstable-T_1))/del_Hc
    Y_O2_2 = Y_O2_1 - (nu*Cp*(T_unstable-T_1))/del_Hc
    Y_N2_2 = 1 - Y_Fuel_2 - Y_O2_2

    # checking for blowout condition
    if np.isnan(T_cold) or np.isnan(T_unstable) or np.isnan(T_hot):
        print("Blow off limit reached for mass flow rate = ", m, " kg/s")  # printing message for blowout
        break

    # appending values in lists
    T_list_1.append(T_cold)
    T_list_2.append(T_unstable)
    T_list_3.append(T_hot)
    Y_Fuel_list.append(Y_Fuel_2)
    Y_O2_list.append(Y_O2_2)
    Y_N2_list.append(Y_N2_2)

# removing initial guess to avoid abrupt jumps in graphs
T_list_1.remove(T_list_1[0])
T_list_2.remove(T_list_2[0])
T_list_3.remove(T_list_3[0])

mass = mass[0:len(T_list_1)]  # making the length of mass flow rate list same as temperature list

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
