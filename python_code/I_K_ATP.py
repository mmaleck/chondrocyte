import numpy as np
from scipy.integrate.odepack import odeint
import matplotlib.pyplot as plt
import functions
from chondrocyte import Voltage_clamp
from params import params_dict
import matplotlib as mpl

"""
The code is used to create Figure 3 for submitted paper 
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

params_dict.update(Mg_i=10)

# Define initial condition vector
y0 = (params_dict["V_0"], params_dict["Na_i_0"], params_dict["K_i_0"], params_dict["Ca_i_0"], params_dict["H_i_0"], 
      params_dict["Cl_i_0"], params_dict["a_ur_0"], params_dict["i_ur_0"], params_dict["vol_i_0"], 
      params_dict["cal_0"])

# fig = plt.figure()
# ax = fig.add_axes([0, 0, 1, 1])

fig, ax = plt.subplots()

solution1 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current1 = Voltage_clamp(solution1)

ax.plot(VV, current1["I_K_ATP"], label="$\mathrm{[K^{+}]_o}$=%d mM" %(params_dict["K_o"]), color="k")

# plt.plot(t, solution1[:,0], label="Mg_i={}".format(params_dict["Mg_i"]))
params_dict.update(K_o=30)
solution2 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current2 = Voltage_clamp(solution2)

ax.plot(VV, current2["I_K_ATP"], label="$\mathrm{[K^{+}]_o}$=%d mM" %(params_dict["K_o"]), color="b")
# plt.plot(t, solution2[:,0], label="Mg_i={}".format(params_dict["Mg_i"]))

params_dict.update(K_o=70)
solution3 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current3 = Voltage_clamp(solution3)

ax.plot(VV, current3["I_K_ATP"], label="$\mathrm{[K^{+}]_o}$=%d mM" %(params_dict["K_o"]), color="r")
# plt.plot(t, solution3[:,0], label="Mg_i={}".format(params_dict["Mg_i"]))
# plt.title("I_K_ATP current with K_o={}".format(params_dict["K_o"]))

ax.set_xlabel("Membrane Potential [mV]", fontsize=16) 
ax.set_ylabel("$\mathrm{I_{K-ATP}}$ [pA/pF]", fontsize=16)
ax.xaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.yaxis.set_tick_params(which='major', size=14, width=2, direction='out')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.legend(loc='lower right', fontsize=16)
# plt.savefig("Fig3_C.png", bbox_inches='tight')
plt.show()


