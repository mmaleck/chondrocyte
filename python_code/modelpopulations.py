"""
This file is part of the chondrocyte modelling project at Simula
Research Laboratory, Norway. Refer to the files README and LICENSE for
more information about the project as well as terms of distribution.
 
Author : Sofie Fischer, M.M.Malacker
Data created : July, 2021
Python version : 3.8.2
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
from chondrocyte import Voltage_clamp
import functions
from params import params_dict

count = 0
num_model = 20

# set time array
t_final = params_dict["t_final"]
dt = params_dict["dt"]
t = np.linspace(0, t_final, int(t_final/dt))

##sample parameter
#inital ion concentrations
sample_Na_i_0 = np.random.lognormal(mean=np.log(params_dict["Na_i_0"]), sigma=0.15, size=num_model)
sample_K_i_0 = np.random.lognormal(mean=np.log(params_dict["K_i_0"]), sigma=0.15, size=num_model)
sample_Ca_i_0 = np.random.lognormal(mean=np.log(params_dict["Ca_i_0"]), sigma=0.15, size=num_model)

#internal volume
sample_vol_i_0 = np.random.lognormal(mean=np.log(params_dict["vol_i_0"]), sigma=0.15, size=num_model)

#Background leakage conductance
sample_g_Na_b_bar = np.random.lognormal(mean=np.log(params_dict["g_Na_b_bar"]), sigma=0.15, size=num_model)
sample_g_K_b_bar = np.random.lognormal(mean=np.log(params_dict["g_K_b_bar"]), sigma=0.15, size=num_model)

#Constants related to the ultra-rapidly rectifying potassium channel
sample_g_K_ur = np.random.lognormal(mean=np.log(params_dict["g_K_ur"]), sigma=0.15, size=num_model)

#Constants related to the two-pore potassium channel
sample_P_K = np.random.lognormal(mean=np.log(params_dict["P_K"]), sigma=0.15, size=num_model)

#Constants related to the calcium-activated potassium channel
sample_gBK = np.random.lognormal(mean=np.log(params_dict["gBK"]), sigma=0.15, size=num_model)

#Temperature-dependance
#sample_I_NaK_scale = np.random.uniform(low=params_dict["g_K_ur"], high=params_dict["g_K_ur"]*4.65, size=num_model)
figs = []
axs = []
noPlots = 3

for i in range(noPlots):
    fig, ax = plt.subplots()
    figs.append(fig)
    axs.append(ax)

while count < num_model:
    
    params_dict["Na_i_0"] = sample_Na_i_0[count]
    params_dict["K_i_0"] = sample_K_i_0[count]
    params_dict["Ca_i_0"] = sample_Ca_i_0[count]
    params_dict["vol_i_0"] = sample_vol_i_0[count]
    params_dict["g_Na_b_bar"] = sample_g_Na_b_bar[count]
    params_dict["g_K_b_bar"] = sample_g_K_b_bar[count]
    params_dict["g_K_ur"] = sample_g_K_ur[count]
    params_dict["gBK"] = sample_gBK[count]
    #params_dict["I_NaK_scale"] = sample_I_NaK_scale[count]
  
    V_0 = params_dict["V_0"]; Na_i_0 = params_dict["Na_i_0"]; K_i_0 = params_dict["K_i_0"]; Ca_i_0 = params_dict["Ca_i_0"]; H_i_0 = params_dict["H_i_0"]; Cl_i_0 = params_dict["Cl_i_0"]
    a_ur_0 = params_dict["a_ur_0"]; i_ur_0 = params_dict["i_ur_0"]; vol_i_0 = params_dict["vol_i_0"]; cal_0 = params_dict["cal_0"]
    
    y0 = (V_0, Na_i_0, K_i_0, Ca_i_0, H_i_0, Cl_i_0, a_ur_0, i_ur_0, vol_i_0, cal_0)
    #params_dict["I_NaK_bar"] = params_dict["I_NaK_scale"]*70.8253*params_dict["C_m"]/params_dict["C_myo"]
    solution = odeint(functions.rhs, y0, t, args=(params_dict,))

    VV, current_dict = Voltage_clamp(solution)
    
    axs[0].plot(VV,current_dict["I_NaK"])
    axs[1].plot(VV,current_dict["I_NaCa"])
    axs[2].plot(VV, current_dict["I_K_2pore"])
    
    count += 1

plt.show()
    
    
    