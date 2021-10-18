import numpy as np
from scipy.integrate.odepack import odeint
import matplotlib.pyplot as plt
import functions
from chondrocyte import Voltage_clamp
from params import params_dict
import matplotlib as mpl

"""
The code is used to create Figure 4A for submitted paper 
"Probing the putative role of KATP channels and biological variability in a mathematical model of chondrocyte electrophysiology”
"""

mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2

# define time span
params_dict.update(t_final=1000)
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
VV, current1 = Voltage_clamp(solution1)

ax.plot(VV, current1["I_K_ATP"], label="$\mathrm{[Mg^{2+}]_i}$=0.1 mM", color="k")

params_dict.update(Mg_i=1.0)
solution2 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current2 = Voltage_clamp(solution2)

ax.plot(VV, current2["I_K_ATP"], label="$\mathrm{[Mg^{2+}]_i}$=1.0 mM", color="b")

params_dict.update(Mg_i=10)
solution3 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current3 = Voltage_clamp(solution3)

ax.plot(VV, current3["I_K_ATP"], label="$\mathrm{[Mg^{2+}]_i}$=10 mM", color="r")

ax.set_xlabel("Membrane Potential [mV]", fontsize=16) 
ax.set_ylabel("$\mathrm{I_{K-ATP}}$ [pA/pF]", fontsize=16)
ax.xaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.yaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(loc="lower right")
plt.savefig("/Users/keiyamamoto/Desktop/Figure_kei/figure4/Fig4_A.png", bbox_inches='tight')
# plt.show()


