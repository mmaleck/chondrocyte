import numpy as np
from scipy.integrate.odepack import odeint
import matplotlib.pyplot as plt
import functions
from chondrocyte import Voltage_clamp
from params import params_dict
import matplotlib as mpl

"""
The code is used to create Figure 2 for submitted paper 
"Probing the putative role of KATP channels and biological variability in a mathematical model of chondrocyte electrophysiology”
"""


mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2

# define time span
params_dict["t_final"] = 50
t_final = params_dict["t_final"]
params_dict["dt"] = 1e-2
dt = params_dict["dt"]
t = np.linspace(0, t_final, int(t_final/dt))

# Define initial condition vector
y0 = (params_dict["V_0"], params_dict["Na_i_0"], params_dict["K_i_0"], params_dict["Ca_i_0"], params_dict["H_i_0"], 
      params_dict["Cl_i_0"], params_dict["a_ur_0"], params_dict["i_ur_0"], params_dict["vol_i_0"], 
      params_dict["cal_0"])

#set simulation specific parameters
params_dict["I_NaK_scale"] = 1.625
params_dict["I_NaK_bar"] = params_dict["I_NaK_scale"]*70.8253*params_dict["C_m"]/params_dict["C_myo"]

params_dict.update(K_o_0=9.5, Na_o=295)

#solve the ODE system which is imported with chondrocyte
solution23 = odeint(functions.rhs, y0, t, args=(params_dict,))

VV, current23 = Voltage_clamp(solution23)


figs = []
axs = []

for i in range(2):
    fig, ax = plt.subplots()
    figs.append(fig)
    axs.append(ax)

I_Nab23 = np.loadtxt('temp23/IV_INab.txt')
I_NaK23 = np.loadtxt('temp23/IV_INaK.txt')
I_K2pore23 = np.loadtxt('temp23/IV_K2pore.txt')

matlab_I_Nab23 = I_Nab23[:, 1]
matlab_I_NaK23 = I_NaK23[:, 1]
matlab_I_K2pore23 = I_K2pore23[:, 1]

axs[0].plot(VV[500:], current23["I_NaK"][500:], label="$\mathrm{I_{Na,K}}$ (Python)", color="k")
axs[0].plot(VV[500:], matlab_I_NaK23[500:], "g--", label="$\mathrm{I_{Na,K}}$ (MATLAB)")
axs[0].plot(VV[500:], current23["I_K_2pore"][500:], label="$\mathrm{I_{K,2pore}}$ (Python)", color="b")
axs[0].plot(VV[500:], matlab_I_K2pore23[500:], "c--", label="$\mathrm{I_{K,2pore}}$ (MATLAB)")
axs[0].plot(VV[500:], current23["I_Na_b"][500:], label="$\mathrm{I_{Na,b}}$ (Python)", color="r")
axs[0].plot(VV[500:], matlab_I_Nab23[500:], "m--", label="$\mathrm{I_{Na,b}}$ (MATLAB)")
axs[0].set_xlabel("Membrane Potential [mV]", fontsize=16) 
axs[0].set_ylabel("Current density [pA/pF]", fontsize=16)
axs[0].xaxis.set_tick_params(which='major', size=14, width=2, direction='out')
axs[0].yaxis.set_tick_params(which='major', size=14, width=2, direction='out')
axs[0].spines['right'].set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[0].legend(loc='lower right', fontsize=12) 
axs[0].set_title("23$^\circ$C ", fontsize=16)

# plt.savefig("Figure2_A.png",bbox_inches='tight')

params_dict["I_NaK_scale"] = 1.625*4.65
params_dict["I_NaK_bar"] = params_dict["I_NaK_scale"]*70.8253*params_dict["C_m"]/params_dict["C_myo"]

solution37 = odeint(functions.rhs, y0, t, args=(params_dict,))
VV, current37 = Voltage_clamp(solution37)


I_Nab37 = np.loadtxt('temp37/IV_INab.txt')
I_NaK37 = np.loadtxt('temp37/IV_INaK.txt')
I_K2pore37 = np.loadtxt('temp37/IV_K2pore.txt')

matlab_I_Nab37 = I_Nab37[:, 1]
matlab_I_NaK37 = I_NaK37[:, 1]
matlab_I_K2pore37 = I_K2pore37[:, 1]

axs[1].plot(VV[500:], current37["I_NaK"][500:], label="$\mathrm{I_{Na,K}}$ (Python)", color="k")
axs[1].plot(VV[500:], matlab_I_NaK37[500:], "g--", label="$\mathrm{I_{Na,K}}$ (MATLAB)")
axs[1].plot(VV[500:], current37["I_K_2pore"][500:], label="$\mathrm{I_{K,2pore}}$ (Python)", color="b")
axs[1].plot(VV[500:], matlab_I_K2pore37[500:], "c--", label="$\mathrm{I_{K,2pore}}$ (MATLAB)")
axs[1].plot(VV[500:], current37["I_Na_b"][500:], label="$\mathrm{I_{Na,b}}$ (Python)", color="r")
axs[1].plot(VV[500:], matlab_I_Nab37[500:], "m--", label="$\mathrm{I_{Na,b}}$ (MATLAB)")
axs[1].set_xlabel("Membrane Potential [mV]", fontsize=16) 
axs[1].set_ylabel("Current density [pA/pF]", fontsize=16)
axs[1].xaxis.set_tick_params(which='major', size=14, width=2, direction='out')
axs[1].yaxis.set_tick_params(which='major', size=14, width=2, direction='out')
axs[1].spines['right'].set_visible(False)
axs[1].spines['top'].set_visible(False)
# axs[1].legend(loc='lower right', fontsize=10)
axs[1].set_title("37$^\circ$C ", fontsize=16)

# plt.savefig("Figure2_B.png",bbox_inches='tight')
plt.show()