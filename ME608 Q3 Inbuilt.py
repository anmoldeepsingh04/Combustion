"""
This code is an alternate solution used to solve the third question of Fundamentals of Combustion (ME608) group project.
Developed by:
Anmoldeep Singh 180030002
"""
# importing relevant libraries
import numpy as np
from scipy.integrate import odeint
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
Ru = 8.315  # in kJ/kmol-K

# Defining initial conditions
P_0 = P_ini * (compression_ratio ** gamma)  # initial pressure
T_0 = T_ini * (compression_ratio ** (gamma - 1))  # initial temperature
concentration_Oxidiser_0 = ((air_to_fuel/phi)/((air_to_fuel/phi)+1)) * (P_0 / (Ru * T_0))  # initial oxidiser concentration
concentration_Fuel_0 = (P_0 / (Ru * T_0)) - concentration_Oxidiser_0  # initial fuel concentration
concentration_Product_0 = 0  # initial product concentration


# to calculate value of the derivatives
def dx_dt(x, t, *args):

    # assigning values to arguments
    concentration_Fuel, concentration_Oxidiser, concentration_Product, T, P = x

    # calculating the concentrations
    concentration_O2 = 0.21*concentration_Oxidiser
    d_concentration_Fuel_dt = (-6.19*(10**9))*np.exp(-15098/T)*(concentration_Fuel**0.1)*(concentration_O2**1.65)
    d_concentration_Oxidiser_dt = 16*d_concentration_Fuel_dt
    d_concentration_Product_dt = -17*d_concentration_Fuel_dt
    dT_dt = ((-h_Fuel*Ru*T)/((Cp-Ru)*P))*d_concentration_Fuel_dt
    dP_dt = (P_0/T_0)*dT_dt
    return d_concentration_Fuel_dt, d_concentration_Oxidiser_dt, d_concentration_Product_dt, dT_dt, dP_dt


# to calculate value of the dp_dt
def delP_dt(arr1, arr2):
    for i in range(0, len(arr2)-1):
        arr1.append((P_0*(arr2[i+1] - arr2[i]))/(T_0*del_t))
    return arr1


# declaring lists to store values
dP_list = [0]
t_start = 0  # starting of time
t_end = 0.15  # ending of time
time_steps = 10000  # number of steps
del_t = (t_end-t_start)/time_steps  # time step
t = np.linspace(t_start, t_end, time_steps)  # to store the time step values
x0 = concentration_Fuel_0, concentration_Oxidiser_0, concentration_Product_0, T_0, P_0  # initial arguments
sol = odeint(dx_dt, x0, t)  # calling function to solve for quantities

# Plotting the data

# variation of Fuel Concentration with Time
plt.plot(t, sol[:, 0], 'r', lw=2.5, label='Fuel')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of fuel (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Fuel Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Oxidiser Concentration with Time
plt.plot(t, sol[:, 1], 'g', lw=2.5, label='Oxidiser')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of oxidiser (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Oxidiser Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Product Concentration with Time
plt.plot(t, sol[:, 2], 'b', lw=2.5, label='Product')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Concentration of product (kmol/m3)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Product Concentration with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Temperature with Time
plt.plot(t, sol[:, 3], 'r', lw=2.5, label='Temperature')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Temperature (K)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Temperature with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# variation of Pressure with Time
plt.plot(t, sol[:, 4], 'm', lw=2.5, label='Pressure')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('Pressure (atm)', fontsize=20)
plt.yticks(fontsize=15)
plt.title('Variation of Pressure with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

# calculating delP_dt
delP_list = delP_dt(dP_list, sol[:, 3])

# variation of dP_dt with Time
plt.plot(t, delP_list, 'c', lw=2.5, label='dP_dt')
plt.xlabel('Time (s)', fontsize=20)
plt.xticks(fontsize=15)
plt.ylabel('dP_dt (atm/s)', fontsize=20)
plt.yticks(fontsize=15)
plt.yscale('log')
plt.title('Variation of dP_dt with Time', fontsize=20)
plt.legend(fontsize=15)
plt.grid()
plt.show()

print("The final fuel concentration is ", sol[-1, 0], "kmol/m3")
print("The final oxidiser concentration is ", sol[-1, 1], "kmol/m3")
print("The final product concentration is ", sol[-1, 2], "kmol/m3")
print("The final temperature is ", sol[-1, 3], "K")
print("The final pressure is ", sol[-1, 4], "atm")