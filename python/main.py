import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import params
from chondrycyte import rhs

'''
following is just the idea how to solve ode. 
'''

# Define time span
t = (0, 5)

# Define initial condition vector
y0 = (params.V_0, params.Na_i_0, params.K_i_0, params.Ca_i_0, params.H_i_0, params.Cl_i_0, params.a_ur_0, params.i_ur_0, params.vol_i_0, params.cal_0)

# Define parameter vector
g_K_bar = params.g_K_bar
P_K = params.P_K
Gmax = params.Gmax
params = (g_K_bar, P_K, Gmax)

# Call the ODE solver
solution = solve_ivp(rhs, t, y0, args=params, max_step=0.01)

t = solution.t
y = solution.y

# Split up into individual states
Mb, O2, MbO2 = y

plt.plot(t, Mb, label='Myoglobin')
plt.plot(t, O2, label='Oxygen')
plt.plot(t, MbO2, label='Oxymyoglobin')

plt.xlabel('Time')
plt.ylabel('Concentrations')
plt.legend()
plt.show()
