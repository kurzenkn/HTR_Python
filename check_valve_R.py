import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from HTR_fluid_prop_eqs import EQUATIONS

L = 15*25.4/1000    #   m
r = 1588*10**-6     #   m
T = 24.5            #   deg C
eta = EQUATIONS["water_props"]["eta_func"](T)     #   Pa*s
A = np.pi*r**2      #   m^2

R_tube = 2*(8*eta*L/(A*r**2))   #   Pa*s/m^3

dP = np.array([10000, 15000, 20000, 10000, 15000, 20000, 10000, 15000, 20000])
t_array = np.array([27.8, 18, 13.4, 27.12, 19.05, 13.8, 27.9, 18.96,13.9])
Q = 8*(10**-6)/t_array

regress_eq = linregress(dP, Q)
slope = regress_eq.slope        #   m^3/Pa*s
R_h_tot = 1/slope               #   Pa*s/m^3

R_h = (R_h_tot-R_tube)/10**9    #   GPa*s/m^3

print("R_h = ", R_h)

x = np.array([10000, 15000, 20000])
y_pred = regress_eq.slope * x + regress_eq.intercept

plt.scatter(dP, Q, label="Measured data")
plt.plot(x, y_pred, color="red", label="Linear regression")

plt.xlabel("dP (Pa)")
plt.ylabel("Q (m^3/s)")
plt.legend()
plt.grid()
plt.show()
