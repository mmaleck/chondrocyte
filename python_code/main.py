import params
from chondrocyte import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
import os

"""
File name : main.py 
Author : Kei Yamamoto
email : keiya@math.uio.no
Data created : July, 2021
Data last modified : 
Python version : 3.8.2
copyright :
credits : 
license : 
"""

path = os.getcwd()
path = os.path.join(path, "result")

if not os.path.isdir(path):
    print("Create a result folder")
    os.makedirs(path)

def main():
    # Define time span
    t_final = params.t_final
    dt = params.dt
    t = np.linspace(0, t_final, int(t_final/dt))

    # Define initial condition vector
    V_0 = params.V_0; Na_i_0 = params.Na_i_0; K_i_0 = params.K_i_0; Ca_i_0 = params.Ca_i_0; H_i_0 = params.H_i_0; Cl_i_0 = params.Cl_i_0
    a_ur_0 = params.a_ur_0; i_ur_0 = params.i_ur_0; vol_i_0 = params.vol_i_0; cal_0 = params.cal_0

    y0 = (V_0, Na_i_0, K_i_0, Ca_i_0, H_i_0, Cl_i_0, a_ur_0, i_ur_0, vol_i_0, cal_0)

    # Define parameter vector
    g_K_b_bar = params.g_K_b_bar
    P_K = params.P_K
    Gmax = params.Gmax

    parameters = (g_K_b_bar, P_K, Gmax)

    # Call the ODE solver
    solution = odeint(rhs, y0, t, args=parameters)

    # Split up into individual states
    # TODO : Cl_i has slightly higher values than MATLAB output  (by Kei)
    V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal  = np.hsplit(solution, 10)

    ramp_Vm = params.ramp_Vm
    if (ramp_Vm == True):
        V_0 = params.V_0; V_final = params.V_final
        V = V_0 + (V_final - V_0)*t/t_final
     
    longth = V.shape[0]

    # get steady-state ion concentrations
    # TODO : only Ca_i_ss is used ? (by Kei)
    Ca_i_ss = Ca_i[longth-1]
    K_i_ss = K_i[longth-1]
    Na_i_ss = Na_i[longth-1]
    V_RMP_ss = V[longth-1]

    # prepare voltage as an array
    V_step_size = 2501
    VV = np.linspace(-150, 100, V_step_size)

    # create dictionary for saving currents
    current_dict = {"alpha_K_DR" : np.zeros(V_step_size), 
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

    # read/declare some parameters outside for loop 
    C_m = params.C_m; K_o = params.K_o
    Na_i_clamp = params.Na_i_clamp; Q = params.Q 
    E_K = -94.02 # From Zhou, et al
    E_Na = 55.0; E_K = -83 
    enable_I_K_DR = params.enable_I_K_DR; enable_I_NaK = params.enable_I_NaK
    enable_I_NaK = params.enable_I_NaK; enable_I_Ca_ATP = params.enable_I_Ca_ATP
    enable_I_K_2pore = params.enable_I_K_2pore; enable_I_Na_b = params.enable_I_Na_b
    enable_I_K_b = params.enable_I_K_b; enable_I_Cl_b = params.enable_I_Cl_b
    enable_I_leak = params.enable_I_leak; enable_I_K_Ca_act = params.enable_I_K_Ca_act

    for i in range(V_step_size):

        # I_K_DR (printed in pA/pF)
        current_dict["I_K_DR"][i], current_dict["alpha_K_DR"][i] = DelayedRectifierPotassium(VV[i], enable_I_K_DR)
        current_dict["I_K_DR"][i] = (current_dict["I_K_DR"][i])/C_m
        
        # I_Na_K (pA; printed IV pA/pF)
        current_dict["I_NaK"][i]= sodiumPotassiumPump(VV[i], K_o, Na_i_clamp,enable_I_NaK)/C_m

        # I_NaCa (pA; printed IV pA/pF)
        current_dict["I_NaCa"][i] = sodiumCalciumExchanger(VV[i], Ca_i_0, Na_i_clamp, enable_I_NaK)/C_m

        # I_Ca_ATP (pA)
        current_dict["I_Ca_ATP"][i] = calciumPump(Ca_i_ss, enable_I_Ca_ATP)

        # I_K_ATP (pA?) Zhou/Ferrero, Biophys J, 2009
        # TODO: it is complex number in the beginning of iterations. need to fix (by Kei)
        current_dict["I_K_ATP"][i] = potassiumPump(VV[i], 0, K_o, E_K, True)

        # I_K_2pore (pA; pA/pF in print) 
        # modeled as a simple Boltzmann relationship via GHK, scaled to match isotonic K+ data from Bob Clark
        current_dict["I_K_2pore"][i] = twoPorePotassium(VV[i], K_i_0, K_o, Q, enable_I_K_2pore)/C_m

        # I_Na_b (pA; pA/pF in print)
        current_dict["I_Na_b"][i] = backgroundSodium(VV[i], None, E_Na, enable_I_Na_b)/C_m

        # I_K_b (pA; pA/pF in print)
        current_dict["I_K_b"][i] = backgroundPotassium(VV[i], None, None, g_K_b_bar, E_K, enable_I_K_b)/C_m
        
        # I_Cl_b (pA; pA/pF in print)
        current_dict["I_Cl_b"][i] = backgroundChloride(VV[i], None, enable_I_Cl_b)/C_m
        
        # I_leak (pA); not printed, added to I_bg
        current_dict["I_leak"][i] = backgroundLeak(VV[i], enable_I_leak)

        # I_bg (pA; pA/pF in print)
        current_dict["I_bq"][i] = current_dict["I_Na_b"][i] + current_dict["I_K_b"][i] + current_dict["I_Cl_b"][i] + current_dict["I_leak"][i]

        # I_K_Ca_act (new version) (pA; pA/pF in print), with converted Ca_i units for model
        current_dict["I_BK"][i] = calciumActivatedPotassium(VV[i], Ca_i_ss, enable_I_K_Ca_act)/C_m

        # I TRPV4 (pA; pA/pF in print)
        current_dict["I_TRPV4"][i] = TripCurrent(VV[i], True)/C_m

        # I_RMP (pA; pA/pF in print)
        current_dict["I_RMP"][i] = current_dict["I_bq"][i] + current_dict["I_BK"][i] + current_dict["I_K_DR"][i] \
                                   + current_dict["I_NaCa"][i] + current_dict["I_NaK"][i] + current_dict["I_K_2pore"][i]

        # I_total (pA)
        # TODO: not yet verifiyed (by Kei)
        current_dict["I_total"][i] = current_dict["I_NaK"][i]*C_m + current_dict["I_NaCa"][i]*C_m + current_dict["I_Ca_ATP"][i] + \
                                     current_dict["I_K_DR"][i]*C_m +  current_dict["I_K_2pore"][i]*C_m + current_dict["I_K_ATP"][i] + \
                                     current_dict["I_BK"][i] + current_dict["I_Na_b"][i]*C_m + current_dict["I_K_b"][i]*C_m + \
                                     current_dict["I_Cl_b"][i]*C_m + current_dict["I_leak"][i] + current_dict["I_TRPV4"][i]*C_m
    
    slope_G = (current_dict["I_bq"][-1]-current_dict["I_bq"][0])/(VV[-1]-VV[0]) # pA/mV = nS
    R = 1/slope_G # = GOhms
    print("slope_G = {}, R={}".format(slope_G, R))
    
    # this is added for debugging. Easier to debug with ipython than using print statement. feel free to remove.
    from IPython import embed; embed(); exit(1)
    
if __name__ == "__main__":
    main()
    
