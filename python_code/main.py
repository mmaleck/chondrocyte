import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import params
from chondrocyte import *

# Define time span
t_final = params.t_final
dt = params.dt
t = np.linspace(0, t_final, int(t_final/dt))

# Define initial condition vector
y0 = (params.V_0, params.Na_i_0, params.K_i_0, params.Ca_i_0, params.H_i_0, params.Cl_i_0, params.a_ur_0, params.i_ur_0, params.vol_i_0, params.cal_0)

# Define parameter vector
g_K_b_bar = params.g_K_b_bar
P_K = params.P_K
Gmax = params.Gmax
parameters = (g_K_b_bar, P_K, Gmax)

# Call the ODE solver
solution = solve_ivp(rhs, t, y0, args=parameters, max_step=0.01)

t = solution.t
y = solution.y

# Split up into individual states
V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal  = y


