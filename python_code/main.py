import params
from chondrocyte import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint
import os

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
    # solution = solve_ivp(fun=rhs, t_span=(0, int(t_final)), y0=y0, method='RK45', t_eval=t, args=parameters)
    solution = odeint(rhs, y0, t, args=parameters)

    # Split up into individual states
    # Cl_i has slightly higher values than MATLAB output 
    V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal  = np.hsplit(solution, 10)

    ramp_Vm = params.ramp_Vm
    if (ramp_Vm == True):
        V_0 = params.V_0; V_final = params.V_final
        V = V_0 + (V_final - V_0)*t/t_final
        #  elseif 
        #     (clamp_Vm == true)
        #     global V_0, global V_final, global t_final;
        #     for

    longth = V.shape[0]


    # get steady-state ion concentrations
    Ca_i_ss = Ca_i[longth-1]
    K_i_ss = K_i[longth-1]
    Na_i_ss = Na_i[longth-1]
    # %Na_i_ss = 12.0
    V_RMP_ss = V[longth-1]



    V_step_size = 2501
    VV = np.linspace(-150, 100, V_step_size)

    
    # create dictionary for saving currents
    current_dict = {"alpha_K_DR" : np.zeros(V_step_size), 
                    "I_K_DR" : np.zeros(V_step_size), 
                    "I_NaK" : np.zeros(V_step_size),
                    "I_NaCa" : np.zeros(V_step_size), 
                    "I_Ca_ATP" : np.zeros(V_step_size),
                    "I_K_ATP" : np.zeros(V_step_size),

                    }

    C_m = params.C_m
    K_o = params.K_o
    Na_i_clamp = params.Na_i_clamp

    for i in range(V_step_size):

        # I_K_DR (printed in pA/pF)
        current_dict["I_K_DR"][i], current_dict["alpha_K_DR"][i] = DelayedRectifierPotassium(VV[i])
        current_dict["I_K_DR"][i] = (current_dict["I_K_DR"][i])/C_m
        
        # I_Na_K (in pA; printed IV pA/pF)
        current_dict["I_NaK"][i]= sodiumPotassiumPump(VV[i], K_o, Na_i_clamp)/C_m

        # I_NaCa (in pA; printed IV pA/pF)
        current_dict["I_NaCa"][i] = sodiumCalciumExchanger(VV[i], Ca_i_0, Na_i_clamp)/C_m

        # I_Ca_ATP (pA)
        current_dict["I_Ca_ATP"][i] = calciumPump(Ca_i_ss)

        # % I_K_ATP (pA?) Zhou/Ferrero, Biophys J, 2009
        # E_K = -94.02; # From Zhou, et al
        # current_dict["I_K_ATP"][i] = potassiumPump(VV[i], None, K_o, E_K)

    from IPython import embed; embed(); exit(1)
        

if __name__ == "__main__":
    main()
    

 

# RFE : possibly save all the values inside the list and write the list into the txt file at the end.
# ofile8 = os.path.join(path,'IV_K2pore.txt')
# fid8 = open(ofile8,'w')
# ofile9 = os.path.join(path,'voltage.txt')
# fid9 = open(ofile9,'w')
# ofile10 = os.path.join(path,'IV_KATP.txt')
# fid10 = open(ofile10,'w')
# ofile11 = os.path.join(path,'IV_INab.txt')
# fid11 = open(ofile11,'w')
# ofile12 = os.path.join(path,'IV_IKb.txt')
# fid12 = open(ofile12,'w')
# ofile13 = os.path.join(path,'IV_IClb.txt')
# fid13 = open(ofile13,'w')
# ofile14 = os.path.join(path,'IV_INaK.txt')
# fid14 = open(ofile14,'w')
# ofile15 = os.path.join(path,'IV_INaCa.txt')
# fid15 = open(ofile15,'w')
# ofile16 = os.path.join(path,'IV_ICaP.txt')
# fid16 = open(ofile16,'w')
# ofile17 = os.path.join(path,'IV_I_TRPV4.txt')
# fid17 = open(ofile17,'w')
# ofile18 = os.path.join(path,'IV_I_bg.txt')
# fid18 = open(ofile18,'w')
# ofile19 = os.path.join(path,'IV_I_total.txt')
# fid19 = open(ofile19,'w')
# ofile20 = os.path.join(path,'IV_I_K_Ca.txt')
# fid20 = open(ofile20,'w')
# ofile21 = os.path.join(path,'IV_I_K_DR.txt')
# fid21 = open(ofile21,'w')
# ofile22 = os.path.join(path,'IV_act_DR.txt')
# fid22 = open(ofile22,'w')
# ofile23 = os.path.join(path,'RMP_vs_KDR.txt')
# fid23 = open(ofile23,'w')
# ofile24 = os.path.join(path,'IV_TRPV4.txt')
# fid24 = open(ofile24,'w')
# ofile25 = os.path.join(path,'IV_RMP.txt')
# fid25 = open(ofile25,'w')

# fid8.close()
# fid9.close()
# fid10.close()
# fid11.close()
# fid12.close()
# fid13.close()
# fid14.close()
# fid15.close()
# fid16.close()
# fid17.close()
# fid18.close()
# fid19.close()
# fid20.close()
# fid21.close()
# fid22.close()
# fid23.close()
# fid24.close()
# fid25 .close()





