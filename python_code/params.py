from math import sqrt

"""
This file is part of the chondrocyte modelling project at Simula
Research Laboratory, Norway. Refer to the files README and COPYING for
more information about the project as well as terms of distribution.
 
Author : Kei Yamamoto, Sofie Fischer, M.M.Malacker
email : keiya@math.uio.no
Data created : July, 2021
Python version : 3.8.2
"""

params_dict = dict(
  clamp_conc = False,             # Toggle clamping of the internal concentrations (for debugging)
  apply_Vm = False,               # Toggle setting of the membrane voltage
  clamp_Vm = True,
  step_Vm = False,
  ramp_Vm = False,
  clamp_Na_i = True,              # Clamp intracellular concentrations to simulate experimental conditions
  clamp_K_i = True,               # Clamp intracellular concentrations to simulate experimental conditions
  V_final = 100.0,
  # Time-stepping information
  t_final = 50000.0,              # Final time (s)
  dt = 1e0,
  # External concentrations
  Na_o = 295.0,                   # Clamped external sodium concentration (mM/l)
  K_o_0  = 9.5,
  K_o = 9.5,                      # must be same as K_o_0
  step_K_o = False,
  Ca_o = 13,                      # Clamped external calcium concentration (mM/l), 1.8-2.0 for cardiac; 1.5 for synovial fluid, 6-15 for chondrocyte matrix
  H_o    = 10**(-7.4),            # Clamped external hydrogen concentration (mM/l)
  Cl_o   = 140.0,                 # Clamped external chloride concentration (mM/l)
  # Initial conditions
  V_0 = -66.725,                  # Initial membrane potential (mV)
  Na_i_0  =  25.0,                # Initial internal sodium concentration (mM/l), 8.5 for cardiac; 40??? for chondrocyte - using either 12.0 or 20.0
  Na_i_clamp = 25.0,
  K_i_0 = 180.0,                  # Initial internal potassium concentration (mM/l), 5.4 for cardiac; 5 for synovial fluid, 7-12 for chondrocyte matrix      
  Ca_i_0 = 0.00001,               # Initial internal calcium concentration (mM/l), 0.000067 mM/l for cardiac, 0.00001 - 0.00005 for chondrocyte
  H_i_0   = 3.47426156721507e-10, # Initial internal hydrogen concentration (mM/l)
  Cl_i_0  = 60.0,                 # Initial internal chloride concentration (mM/l)"""
  # Initial value of Kv activation
  a_ur_0  = 1.95403736678201e-04, 
  i_ur_0  = 9.99896050539933e-01, 
  cal_0   = 5.08345961310772e-05, # From Nygren, et al 1998
  # Universal constants
  R = 8314.472,                   # Universal gas constant (mJ K^-1 mol^-1)
  T = 310.15,                     # Normal body temperature (K)
  F = 96485.34,                   # Faradays constant (C mol^-1)
  # Charges on each of the ions
  z_Na = 1,                       # Charge on the sodium ion
  z_K  = 1,                       # Charge on the potassium ion
  z_Ca = 2,                       # Charge on the calcium ion
  z_H  = 1,                       # Charge on the calcium ion
  z_Cl  = 1,                      # Charge on the chloride ion
  # Cell parameters
  C_m = 6.3,                      # Membrane capacitance, (pF)
  vol_i_0 = 0.005884,             # Internal volume
  C_myo = 50.0,                   # (pF); myocyte capacitance from Nygren, et al, 1998, for scaling
  # Constants related to external stimulation
  t_cycle = 5.0,                  # Total cycle time (s)
  t_stim = 1.0,                   # Stimulation time/cycle (s)
  I_stim_bar = 0.0,               # Stimulation current magnitude (pA)
  # Background conductances
  g_Na_b_bar = 0.1,               # Background sodium leakage conductance (pS)
  I_Na_b_scale = 1.0,           
  g_K_b_bar = 0.07,               # Background potassium leakage conductance (pS), No background K in Nygren1998
  g_Cl_b_bar = 0.05,              # Background chloride leakage conductance (pS)
  g_leak = 0.0,                   # Background seal leakage - set to 0.5
  I_NaK_scale = 1.625,
  K_NaK_K = 2.1,                  # (mmol/L) Nygren, et al, 1998
  K_NaK_Na = 17.5,                # (mmol/L) Nygren, et al, 1998
  K_NaCa = 0.0374842,             # pA/(mmol/L)4 (Nygren1998) - 0.04 - default
  gamma_Na = 0.45,                # same, dimensionless
  d_NaCa = 0.0003,                # same, (mmol/L)^4
  # Constants related to the sodium-hydrogen exchanger
  n_H = 1,
  m_H = 3,
  K_H_i_mod = 3.07e-5,
  K_H_o_mod = 4.8e-7,
  k1_p = 10.5,
  k1_m = 0.201,
  k2_p = 15.8,
  k2_m = 183,
  K_Na_i = 16.2,
  K_Na_o = 195,
  K_H_i = 6.05e-4,
  K_H_o = 1.62e-3,
  N_NaH_channel = 4899,
  I_NaH_scale = 1.00,
  # Constants related to the calcium pump
  I_Ca_ATP_bar = 6.0,             # ICaP = 4.0 (pA), Nygren1998
  I_Ca_ATP_scale = 1.0, 
  k_Ca_ATP = 0.0002,              # kCaP = 0.0002 (mmol/L)
  # Constants related to the ultra-rapidly rectifying potassium channel
  g_K_ur = 1.4623,
  # Constants related to I_K_DR
  i_K_DR = 1.0,
  act_DR_shift = 10.0,
  # Constants related to the two-pore potassium channel
  P_K = 3.1e-6*sqrt(5/140),
  Q = 0.1,                        # Based on experimental data from Bob; slope of IV in isopotential recording conditions
  I_K_2pore_0 = 0.0,
  I_K_2pore_scale = 1.35, 
  # Constants related to the calcium-activated potassium channel
  Zj = 0.70,
  Vhj = 250,
  ZL = 0.1,
  L0 = 12e-6,
  KDc = 3e-6,
  C = 8,
  D = 25,
  E = 2.4,
  Gmax = 3.8*2.4,
  N_channel = 1.0,
  E_K_Ca_act = 42,
  gBK = 2.50,
  # Constants related to the TRP1 channel
  g_TRP1 = 1.e-4*0.75*10/4,
  a_TRP1 = 80,
  b_TRP1 = -1000,
  # Constants related to the TRPV4 channel
  g_TRPV4 = 0.00046875*0.05, # (nS)
  a_TRPV4 = 80,
  b_TRPV4 = -1000,
  E_Na = 55.0
)

# Constants related to the sodium-potassium pump, pA, from Nygren et al, 1998
params_dict["I_NaK_bar"] = params_dict["I_NaK_scale"]*70.8253*params_dict["C_m"]/params_dict["C_myo"]

# Constants related to the sodium-calcium exchanger
params_dict["NCX_scale"] = params_dict["C_m"] / params_dict["C_myo"]

# Constants related to I_K_DR, given in nS/pF per Bob Clark communication 3/13/16, assuming C_m as above
params_dict["g_K_DR"] = 0.0289*params_dict["C_m"]
