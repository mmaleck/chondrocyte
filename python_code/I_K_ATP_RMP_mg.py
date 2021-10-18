import numpy as np
from scipy.integrate.odepack import odeint
import matplotlib.pyplot as plt
import functions
from chondrocyte import Voltage_clamp
from params import params_dict
import matplotlib as mpl

"""
The code is used to create Figure 4B for submitted paper 
"Probing the putative role of KATP channels and biological variability in a mathematical model of chondrocyte electrophysiology”
"""


mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2


# define time span
params_dict.update(t_final=180)
t_final = params_dict["t_final"]
dt = params_dict["dt"]
t = np.linspace(0, t_final, int(t_final/dt))

params_dict.update(Mg_i=1)

# Define initial condition vector
y0 = (params_dict["V_0"], params_dict["Na_i_0"], params_dict["K_i_0"], params_dict["Ca_i_0"], params_dict["H_i_0"], 
      params_dict["Cl_i_0"], params_dict["a_ur_0"], params_dict["i_ur_0"], params_dict["vol_i_0"], 
      params_dict["cal_0"])

fig, ax = plt.subplots()

params_dict.update(K_o_0=7, Mg_i=0.1)

solution1 = odeint(functions.rhs, y0, t, args=(params_dict,))
ax.plot(t, solution1[:,0], label="$\mathrm{[Mg^{2+}]_i}$=0.1 mM", color="k")

params_dict.update(Mg_i=1.0)
solution2 = odeint(functions.rhs, y0, t, args=(params_dict,))
ax.plot(t, solution2[:,0], label="$\mathrm{[Mg^{2+}]_i}$=1.0 mM", color="b")

params_dict.update(Mg_i=10)
solution3 = odeint(functions.rhs, y0, t, args=(params_dict,))
ax.plot(t, solution3[:,0], label="$\mathrm{[Mg^{2+}]_i}$=10 mM", color="r")
ax.set_xlabel("Time [s]", fontsize=16) 
ax.set_ylabel("Membrane Potential [mV]", fontsize=16)
ax.xaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.yaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(loc='upper right')
# plt.savefig("Fig4_B.png", bbox_inches='tight')
plt.show()


