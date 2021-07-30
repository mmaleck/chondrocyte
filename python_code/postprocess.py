from pickle import load
from os import path
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from params import params_dict
'''
Postprocess file
Exmaple usage : python postprocess.py --case="result/1" --

current_name : "I_K_DR" 
               "I_NaK" 
               "I_NaCa"
               "I_Ca_ATP"
               "I_K_ATP" 
               "I_K_2pore"
               "I_Na_b" 
               "I_K_b" 
               "I_Cl_b" 
               "I_leak" 
               "I_bq" 
               "I_BK" 
               "I_TRPV4" 
               "I_RMP" 
               "I_total"

Author : Kei Yamamoto
email : keiya@math.uio.no
Data created : July, 2021
Python version : 3.8.2
'''

parser = ArgumentParser()
parser.add_argument('--case', type=str, default="result/1", help="Path to simulation results",
                        metavar="PATH")

args = parser.parse_args()

current_file = open(path.join(args.case,"current.pkl"), "rb")
current = load(current_file)
current_file.close()

ode_solution_file = open(path.join(args.case,"ode_solution.pkl"), "rb")
ode_solution = load(ode_solution_file)
ode_solution_file.close()
V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal  = np.hsplit(ode_solution, 10) 


# Example of plotting 
# I_V plot
V_step_size = params_dict["V_step_size"]
VV = np.linspace(params_dict["V_start"], params_dict["V_end"], V_step_size)
plt.plot(VV, current["I_NaK"])
plt.xlabel("voltage(mV)")
plt.ylabel("I,NaK")
plt.legend()
plt.show()

"""
You can add codes from here to create plot and play around with the simulation results
"""
