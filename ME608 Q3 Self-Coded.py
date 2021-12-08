"""
This code is used to solve the third question of Fundamentals of Combustion (ME608) group project.
Developed by:
Anmoldeep Singh 180030002
"""
# importing relevant libraries
import numpy as np
import matplotlib.pyplot as plt

# Defining required constants
Mol_mass_Fuel, Mol_mass_Oxidiser, Mol_mass_Product = 29.0, 29.0, 29.0  # in kg/kmol
Cp = 1200  # in J/kg-K
h_Fuel = 4.0*(10**7)  # in J/kg
air_to_fuel = 16.0  # air to fuel ratio
phi = 1  # equivalence ratio
gamma = 1.4  # specific heat ratio
T_ini = 300  # in K
P_ini = 1  # in atm
compression_ratio = 10
x = 0.1
y = 1.65
Ru = 8.315  # in kJ/kmol-K
del_t = 0.00001  # time step

# Defining initial conditions
P_0 = P_ini * (compression_ratio ** gamma)  # initial pressure
T_0 = T_ini * (compression_ratio ** (gamma - 1))  # initial temperature
concentration_Oxidiser_0 = ((air_to_fuel/phi)/((air_to_fuel/phi)+1)) * (P_0 / (Ru * T_0))  # initial oxidiser concentration
concentration_Fuel_0 = (P_0 / (Ru * T_0)) - concentration_Oxidiser_0  # initial fuel concentration
concentration_Product_0 = 0  # initial product concentration

# creating lists to store values
concentration_Fuel = [concentration_Fuel_0]
concentration_Oxidiser = [concentration_Oxidiser_0]
concentration_Product = [concentration_Product_0]
T = [T_0]
P = [P_0]
del_P = [0]

# initializing values
concentration_Fuel_old = concentration_Fuel_0
concentration_Oxidiser_old = concentration_Oxidiser_0
concentration_Product_old = concentration_Product_0
T_old = T_0
P_old = P_0
t = [0]
time = 0

while concentration_Fuel_old > 10e-10:  # iterating until concentration of fuel is above accepted limit

    # calculating the concentrations
    concentration_O2 = concentration_Oxidiser_old*0.21
    concentration_Fuel_new = concentration_Fuel_old - del_t*(6.19*(10**9))*np.exp(-15098/T_old)*(concentration_Fuel_old**x)*(concentration_O2**y)
    concentration_Oxidiser_new = 16*concentration_Fuel_new + concentration_Oxidiser_0 - 16*concentration_Fuel_0
    concentration_Product_new = -17*concentration_Fuel_new + concentration_Product_0 + 17*concentration_Fuel_0

    # calculating new temperature
    T_new = T_0 + ((Ru*h_Fuel*T_0*concentration_Fuel_0)/((Cp-Ru)*P_0)) - ((Ru*h_Fuel*T_0*concentration_Fuel_new)/((Cp-Ru)*P_0))

    # calculating new pressure
    P_new = P_old + (P_0/T_0)*(T_new-T_old)

    # calculating new del_P
    del_P_new = (P_0*(T_new - T_old))/(T_0*del_t)

    # updating time step
    time += del_t

    # appending values in lists
    concentration_Fuel.append(concentration_Fuel_new)
    concentration_Oxidiser.append(concentration_Oxidiser_new)
    concentration_Product.append(concentration_Product_new)
    T.append(T_new)
    P.append(P_new)
    del_P.append(del_P_new)
    t.append(time)

    # updating values
    concentration_Fuel_old = concentration_Fuel_new
    concentration_Oxidiser_old = concentration_Oxidiser_new
    concentration_Product_old = concentration_Product_new
    T_old = T_new
    P_old = P_new


# Plotting the data

# variation of Fuel Concentration with Time
plt.plot(t, concentration_Fuel, 'r', lw=2.5, label='Fuel')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of fuel (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Fuel Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Oxidiser Concentration with Time
plt.plot(t, concentration_Oxidiser, 'g', lw=2.5, label='Oxidiser')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of oxidiser (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Oxidiser Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Product Concentration with Time
plt.plot(t, concentration_Product, 'b', lw=2.5, label='Product')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of product (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Product Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Temperature with Time
plt.plot(t, T, 'r', lw=2.5, label='Temperature')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Temperature (K)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Temperature with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Pressure with Time
plt.plot(t, P, 'm', lw=2.5, label='Pressure')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Pressure (atm)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Pressure with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of dP_dt with Time
plt.plot(t, del_P, 'c', lw=2.5, label='dP_dt')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('dP_dt (atm/s)', fontsize=20)
plt.yticks(fontsize=15)
plt.yscale('log')
plt.title('Variation of dP_dt with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

print("The final fuel concentration is ", concentration_Fuel[-1], "kmol/m3")
print("The final oxidiser concentration is ", concentration_Oxidiser[-1], "kmol/m3")
print("The final product concentration is ", concentration_Product[-1], "kmol/m3")
print("The final temperature is ", T[-1], "K")
print("The final pressure is ", P[-1], "atm")
