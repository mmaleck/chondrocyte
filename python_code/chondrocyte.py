import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
import os
import pickle
import functions
from params import params_dict

"""
This file is part of the chondrocyte modelling project at Simula
Research Laboratory, Norway. Refer to the files README and COPYING for
more information about the project as well as terms of distribution.
 
Author : Kei Yamamoto, Sofie Fischer
email : keiya@math.uio.no
Data created : July, 2021
Python version : 3.8.2
"""

# Apply voltage clamp
def Voltage_clamp(solution, paramter):

    # prepare voltage as an array
    V_step_size = paramter["V_step"]
    VV = np.linspace(paramter["V_start"], paramter["V_end"], V_step_size)
    
    # get steady-state ion concentrations
    Ca_i_ss = solution[:,3][solution.shape[0]-1]
    Na_i_ss = solution[:,1][solution.shape[0]-1]

    # create dictionary for saving currents
    current_dict = {
                    "I_K_DR" : np.zeros(V_step_size), 
                    "I_NaK" : np.zeros(V_step_size),
                    "I_NaCa" : np.zeros(V_step_size), 
                    "I_Ca_ATP" : np.zeros(V_step_size),
                    "I_K_ATP" : np.zeros(V_step_size),
                    "I_K_2pore": np.zeros(V_step_size),
                    "I_Na_b" : np.zeros(V_step_size),
                    "I_K_b" : np.zeros(V_step_size),
                    "I_Cl_b" : np.zeros(V_step_size),
                    "I_leak" : np.zeros(V_step_size),
                    "I_bq" : np.zeros(V_step_size),
                    "I_BK" : np.zeros(V_step_size),
                    "I_TRPV4" : np.zeros(V_step_size),
                    "I_RMP" : np.zeros(V_step_size), 
                    "I_total" : np.zeros(V_step_size)
                    }

    # declear parameters
    C_m = paramter["C_m"]; K_o=paramter["K_o"]; Na_i_clamp = paramter["Na_i_clamp"]
    Ca_i_0 = paramter["Ca_i_0"]; K_i_0 = paramter["K_i_0"]; Q=paramter["Q"]
    E_Na=paramter["E_Na"]; g_K_b_bar=paramter["g_K_b_bar"]; temp=paramter["temp"]

    for i in range(V_step_size):

        # I_K_DR (printed in pA/pF)
        current_dict["I_K_DR"][i] = functions.DelayedRectifierPotassium(V=VV[i], enable_I_K_DR=True)/C_m

        # I_Na_K (pA; printed IV pA/pF)
        current_dict["I_NaK"][i]= functions.sodiumPotassiumPump(V=VV[i], K_o=K_o, Na_i_0=Na_i_clamp, enable_I_NaK=True)/C_m

        # I_NaCa (pA; printed IV pA/pF)
        current_dict["I_NaCa"][i] = functions.sodiumCalciumExchanger(V=VV[i], Ca_i=Ca_i_0, Na_i_0=Na_i_clamp, enable_I_NaCa=True)/C_m

        # I_Ca_ATP (pA)
        current_dict["I_Ca_ATP"][i] = functions.calciumPump(Ca_i=Ca_i_ss, enable_I_Ca_ATP=True)

        # I_K_ATP (pA?) Zhou/Ferrero, Biophys J, 2009
        current_dict["I_K_ATP"][i] = functions.potassiumPump(V=VV[i], K_i=None, K_o=K_o,E_K=-94.02, Na_i=Na_i_ss, temp=temp, enable_I_K_ATP=True)

        # I_K_2pore (pA; pA/pF in print) 
        # modeled as a simple Boltzmann relationship via GHK, scaled to match isotonic K+ data from Bob Clark
        current_dict["I_K_2pore"][i] = functions.twoPorePotassium(V=VV[i], K_i_0=K_i_0, K_o=K_o, Q=Q, enable_I_K_2pore=True)/C_m

        # I_Na_b (pA; pA/pF in print)
        current_dict["I_Na_b"][i] = functions.backgroundSodium(V=VV[i], Na_i=None, E_Na=E_Na, enable_I_Na_b=True)/C_m

        # I_K_b (pA; pA/pF in print)
        current_dict["I_K_b"][i] = functions.backgroundPotassium(V=VV[i], K_i=None, K_o=None, g_K_b_bar=g_K_b_bar, E_K=-83, enable_I_K_b=True)/C_m
        
        # I_Cl_b (pA; pA/pF in print)
        current_dict["I_Cl_b"][i] = functions.backgroundChloride(V=VV[i], Cl_i=None, enable_I_Cl_b=True)/C_m
        
        # I_leak (pA); not printed, added to I_bg
        current_dict["I_leak"][i] = functions.backgroundLeak(V=VV[i], enable_I_leak=False)

        # I_bg (pA; pA/pF in print)
        current_dict["I_bq"][i] = current_dict["I_Na_b"][i] + current_dict["I_K_b"][i] + current_dict["I_Cl_b"][i] + current_dict["I_leak"][i]

        # I_K_Ca_act (new version) (pA; pA/pF in print), with converted Ca_i units for model
        current_dict["I_BK"][i] = functions.calciumActivatedPotassium(V=VV[i], Ca_i=Ca_i_0, enable_I_K_Ca_act=True)/C_m

        # I TRPV4 (pA; pA/pF in print)
        current_dict["I_TRPV4"][i] = functions.TripCurrent(V=VV[i], enable_I_TRPV4=True)/C_m

        # I_RMP (pA; pA/pF in print)
        current_dict["I_RMP"][i] = current_dict["I_bq"][i] + current_dict["I_BK"][i] + current_dict["I_K_DR"][i] \
                                + current_dict["I_NaCa"][i] + current_dict["I_NaK"][i] + current_dict["I_K_2pore"][i]

        # I_total (pA)
        current_dict["I_total"][i] = current_dict["I_NaK"][i]*C_m + current_dict["I_NaCa"][i]*C_m + current_dict["I_Ca_ATP"][i] + \
                                    current_dict["I_K_DR"][i]*C_m +  current_dict["I_K_2pore"][i]*C_m + current_dict["I_K_ATP"][i] + \
                                    current_dict["I_BK"][i] + current_dict["I_Na_b"][i]*C_m + current_dict["I_K_b"][i]*C_m + \
                                    current_dict["I_Cl_b"][i]*C_m + current_dict["I_leak"][i] + current_dict["I_TRPV4"][i]*C_m 
    
    # slope_G = (current_dict["I_bq"][-1]-current_dict["I_bq"][0])*C_m/(VV[-1]-VV[0]) # pA/mV = nS
    # R = 1/slope_G # = GOhms

    return VV, current_dict

# Following code will create a folder called "result" when you first run the simulation 
# and also create a folder 1 inside "result" folder.
# If there already exists "result", it will create a new folder based on the number of folders that exits inside "result"
path = os.getcwd()
newfolder = os.path.join(path, "result")

if not os.path.exists(newfolder):
    newfolder = os.path.join(newfolder, '1')
    os.makedirs(newfolder)
else:
    previous = [f for f in os.listdir(newfolder) if not f.startswith('.')]
    previous = max(map(eval, previous)) if previous else 0
    newfolder = os.path.join(newfolder, str(previous + 1))
    os.makedirs(newfolder)

# Define time span
t_final = params_dict["t_final"]
dt = params_dict["dt"]
t = np.linspace(0, t_final, int(t_final/dt))

# Define initial condition vector
y0 = (params_dict["V_0"], params_dict["Na_i_0"], params_dict["K_i_0"], params_dict["Ca_i_0"], params_dict["H_i_0"], 
      params_dict["Cl_i_0"], params_dict["a_ur_0"], params_dict["i_ur_0"], params_dict["vol_i_0"], 
      params_dict["cal_0"])

# Call the ODE solver
solution_ode = odeint(functions.rhs, y0, t, args=(params_dict,))

# CaLL Voltage Clamp
VV, current_dict = Voltage_clamp(solution_ode, params_dict)

# save ode_solution
with open(os.path.join(newfolder, 'ode_solution.pkl'), 'wb') as file:
    pickle.dump(solution_ode, file)

# save current
with open(os.path.join(newfolder, 'current.pkl'), 'wb') as file:
    pickle.dump(current_dict, file)

# save parameters as text file
with open(os.path.join(newfolder, 'params.txt'), 'w') as par:
    for key, value in params_dict.items(): 
        par.write('%s: %s\n' % (key, value))
