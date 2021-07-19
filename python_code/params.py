from math import sqrt

'''
  This file is part of the chondrocyte modelling project at Simula
  Research Laboratory, Norway. Refer to the files README and COPYING for
  more information about the project as well as terms of distribution.
 
  Copyright (C) 2010--2016 M. M. Maleckar & Harish Narayanan
  Licensed under the GNU GPL Version 3
  The solution fields and parameters have the following units:
 
  Voltage: mV
  Current: pa
  Time: s
  Concentration: mM/l
  Conductance: pS except where noted -- nS? 3/16/16
  Capacitance: pF
  Current RMP ~-33.25 mV, Na_i ~5mM and K_i ~125 mM
'''

params_dict = dict(
  enable_parest = False,
  clamp_conc = False,
  apply_Vm = False,
  clamp_Vm = True,
  step_Vm = False,
  ramp_Vm = False,
  clamp_Na_i = True,
  clamp_K_i = True,
  V_final = 100.0,
  t_final = 50000.0,       # Final time (s)
  dt = 1e0,
  Na_o = 295.0,
  K_o_0  = 9.5,
  K_o = 9.5, # must be same as K_o_0
  step_K_o = False,
  Ca_o = 13,
  H_o    = 10**(-7.4),
  Cl_o   = 140.0,
  V_0 = -66.725,
  Na_i_0  =  25.0,
  Na_i_clamp = 25.0,
  K_i_0 = 180.0,    #Initial internal potassium concentration (mM/l)
  Ca_i_0 = 0.00001, #Initial internal calcium concentration (mM/l)
  H_i_0   = 3.47426156721507e-10, #Initial internal hydrogen concentration (mM/l)
  Cl_i_0  = 60.0, #Initial internal chloride concentration (mM/l)"""
  a_ur_0  = 1.95403736678201e-04, #Initial value of Kv activation
  i_ur_0  = 9.99896050539933e-01, 
  cal_0   = 5.08345961310772e-05,
  R = 8314.472, # Universal gas constant (mJ K^-1 mol^-1)
  T = 310.15,   # Normal body temperature (K)
  F = 96485.34, # Faradays constant (C mol^-1)
  z_Na = 1,         # Charge on the sodium ion
  z_K  = 1,        # Charge on the potassium ion
  z_Ca = 2,        # Charge on the calcium ion
  z_H  = 1,         # Charge on the calcium ion
  z_Cl  = 1,       # Charge on the chloride ion
  C_m = 6.3,           # Membrane capacitance, (pF)
  vol_i_0 = 0.005884, # Internal volume
  C_myo = 50.0,        # (pF); myocyte capacitance from Nygren, et al, 1998, for scaling
  t_cycle = 5.0,    # Total cycle time (s)
  t_stim = 1.0,     # Stimulation time/cycle (s)
  I_stim_bar = 0.0, # Stimulation current magnitude (pA)
  g_Na_b_bar = 0.1,
  I_Na_b_scale = 1.0, # for temp 23C
  g_K_b_bar = 0.07,
  g_Cl_b_bar = 0.05,
  g_leak = 0.0,
  I_NaK_scale = 1.625,
  K_NaK_K = 2.1,
  K_NaK_Na = 17.5,
  K_NaCa = 0.0374842, #pA/(mmol/L)4 (Nygren1998) - 0.04 - default
  gamma_Na = 0.45, #same, dimensionless
  d_NaCa = 0.0003, #same, (mmol/L)^4
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
  I_NaH_scale = 1.00, # for temp 23C
  I_Ca_ATP_bar = 6.0, #ICaP = 4.0 (pA), Nygren1998
  I_Ca_ATP_scale = 1.0, # for temp 23C
  k_Ca_ATP = 0.0002, #kCaP = 0.0002 (mmol/L)
  g_K_ur = 1.4623,
  i_K_DR = 1.0,
  act_DR_shift = 10.0,
  P_K = 3.1e-6*sqrt(5/140),
  Q = 0.1, #Based on experimental data from Bob; slope of IV in isopotential recording conditions
  I_K_2pore_0 = 0.0,
  I_K_2pore_scale = 1.35, # for temp 23C
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
  g_TRP1 = 1.e-4*0.75*10/4,
  a_TRP1 = 80,
  b_TRP1 = -1000,
  g_TRPV4 = 0.00046875*0.05, # (nS)
  a_TRPV4 = 80,
  b_TRPV4 = -1000,
  E_Na = 55.0
)

params_dict["I_NaK_bar"] = params_dict["I_NaK_scale"]*70.8253*params_dict["C_m"]/params_dict["C_myo"]
params_dict["NCX_scale"] = params_dict["C_m"] / params_dict["C_myo"]
params_dict["g_K_DR"] = 0.0289*params_dict["C_m"]



# """Toggle estimation of parameters"""
# # TODO : not used anywhere ?
# enable_parest = False

# """Toggle clamping of the internal concentrations (for debugging)"""
# clamp_conc = False

# """Toggle setting of the membrane voltage"""
# apply_Vm = False

# """If apply_Vm is true, then one of clamp_Vm, ramp_Vm and step_Vm must be
# true to define what voltage is to be applied"""
# clamp_Vm = True
# step_Vm = False
# ramp_Vm = False
# # V_final = 90.0  # Final value of membrane voltage when ramped (mV) # TODO: not 100 ?
# V_final = 100.0

# """Clamp intracellular concentrations to simulate experimental conditions"""
# clamp_Na_i = True
# clamp_K_i = True

# """=====Toggle individual currents====="""

# """Background"""
# # enable_I_Na_b = True
# # enable_I_K_b = True
# # enable_I_Cl_b = True
# # enable_I_leak = False

# """Pumps and exchangers"""
# # enable_I_NaK = True
# # enable_I_NaCa = True
# # enable_I_NaH = True
# # enable_I_Ca_ATP = True

# """Potassium currents"""
# # enable_I_K_ur = False     #Leave turned off; redundant as this was replaced by I_K_DR, below
# # enable_I_K_2pore = True   #Basically no effect on RMP in current form
# # enable_I_K_Ca_act = True  #Little effect on RMP in current form
# # enable_I_K_ATP = False    #Large effect on RMP when turned on
# # enable_I_K_DR = True      #Large effect on RMP

# """Other currents"""
# # enable_I_TRP1 = False
# # enable_I_TRPV4 = False

# # """Definitely not used"""
# # enable_I_ASIC = False
# # enable_I_TRP2 = False
# # enable_I_stim = False

# """Time-stepping information"""
# t_final = 50000.0       # Final time (s)
# dt = 1e0
# # t_final = 60
# # dt = 1e-2                   

# """===External concentrations==="""

# """ Na_o  : 130 for cardiac; 140 for synovial fluid, 240-250 for chondrocyte matrix    
# Clamped external sodium concentration (mM/l)"""
# # Na_o = 140
# Na_o = 295.0

# """Clamped external potassium concentration (mM/l)
# K_o_0, K_o : 5.4 for cardiac; 5 for synovial fluid, 7-12 for chondrocyte matrix      
# """
# K_o_0  = 9.5
# # K_o_0 = 5.0
# K_o = K_o_0
# step_K_o = False

# """ Ca_o : 1.8-2.0 for cardiac; 1.5 for synovial fluid, 6-15 for chondrocyte matrix      
# Clamped external calcium concentration (mM/l)"""
# # Ca_o   = 1.5
# Ca_o = 13

# """Clamped external hydrogen concentration (mM/l)"""
# H_o    = 10**(-7.4) 

# """Clamped external chloride concentration (mM/l)"""
# Cl_o   = 140.0        

# """=====Initial conditions====="""

# """Initial membrane potential (mV), -100, -6.17150137799178e+01 as option"""
# # V_0 = -40.0
# V_0 = -66.725

# """Initial internal sodium concentration (mM/l)
# Na_i_0 : 8.5 for cardiac; 40??? for chondrocyte - using either 12.0 or 20.0
# Na_i_0 = 1.22582880260390e+00 """
# Na_i_0  =  25.0 
# # Na_i_0 = 12.0
# Na_i_clamp = 25.0

# """Initial internal potassium concentration (mM/l)
# K_i_0 : 120-140 for cardiac; 120-140 for chondrocyte
# K_i_0 = 1.16805063818441e+02"""
# K_i_0 = 180.0
# # K_i_0 = 140 

# """Initial internal calcium concentration (mM/l)
# Ca_i_0 : 0.000067 mM/l for cardiac, 0.00001 - 0.00005 for chondrocyte
# Ca_i_0 = 0.0015
# Ca_i_0 = 1.20992489429946e-07"""
# Ca_i_0 = 0.00001
# # Ca_i_0 = 0.00005

# """Initial internal hydrogen concentration (mM/l)"""
# H_i_0   = 3.47426156721507e-10

# """Initial internal chloride concentration (mM/l)"""
# Cl_i_0  = 60.0 

# """Initial value of Kv activation"""
# a_ur_0  = 1.95403736678201e-04 
# i_ur_0  = 9.99896050539933e-01 
# cal_0   = 5.08345961310772e-05
# cal_0 = 0.0275 #From Nygren, et al 1998

# """Universal constants"""
# R = 8314.472 # Universal gas constant (mJ K^-1 mol^-1)
# T = 310.15   # Normal body temperature (K)
# F = 96485.34 # Faradays constant (C mol^-1)

# """Charges on each of the ions"""
# z_Na = 1         # Charge on the sodium ion
# z_K  = 1         # Charge on the potassium ion
# z_Ca = 2         # Charge on the calcium ion
# z_H  = 1         # Charge on the calcium ion
# z_Cl  = 1        # Charge on the chloride ion

# """Cell parameters"""
# C_m = 6.3           # Membrane capacitance, (pF)
# vol_i_0 = 0.005884  # Internal volume
# C_myo = 50.0        # (pF); myocyte capacitance from Nygren, et al, 1998, for scaling

# """Constants related to external stimulation"""
# t_cycle = 5.0    # Total cycle time (s)
# t_stim = 1.0     # Stimulation time/cycle (s)
# I_stim_bar = 0.0 # Stimulation current magnitude (pA)

# """======Background conductances====="""

# """Background sodium leakage conductance (pS)
# #g_Na_b_bar = 0.30 # gB,Na = 0.060599 nS, Nygren1998"""
# g_Na_b_bar = 0.1
# I_Na_b_scale = 1.0 # for temp 23C
# # I_Na_b_scale = 1.29 # for temp 37C


# """Background potassium leakage conductance (pS), No background K in Nygren1998
# g_K_b_bar = 0.65  # No background K in Nygren1998
# g_K_b_bar = 0.07 + 0.05
# """
# g_K_b_bar = 0.07


# """Background chloride leakage conductance (pS)
# g_Cl_b_bar = 0.40 """
# g_Cl_b_bar = 0.05

# """Background seal leakage - set to 0.5""" # TODO : not set to 0.5 ? (by Kei)
# g_leak = 0.0 

# """Constants related to the sodium-potassium pump, pA""" # TODO : defined again below, should fix or combine (by Kei)
# # I_NaK_bar = 93.8 #82 (pA) 
# # I_NaK_bar = 70.8253*C_m/C_myo #pA, from Nygren et al, 1998, default

# """EXPERIMENT for INaK @ physiological temperature"""
# # I_NaK_scale = 4.65*1.625                  #scaled for q10 = 3 per Wayne's request, for temp 37C
# # I_NaK_scale = 1.00 # for temp 23
# I_NaK_scale = 1.625
# I_NaK_bar = I_NaK_scale*70.8253*C_m/C_myo     #(pA), from Nygren et al, 1998
# # K_NaK_K = 1.0                               #(mmol/L) Nygren, et al, 1998
# # K_NaK_K = 2.2
# K_NaK_K = 2.1
# # K_NaK_Na = 11.0                             #(mmol/L) Nygren, et al, 1998
# # K_NaK_Na = 12.0
# K_NaK_Na = 17.5

# """Constants related to the sodium-calcium exchanger"""
# K_NaCa = 0.0374842 #pA/(mmol/L)4 (Nygren1998) - 0.04 - default
# # K_NaCa = 0.077472
# # K_NaCa = 0.05
# gamma_Na = 0.45 #same, dimensionless
# d_NaCa = 0.0003 #same, (mmol/L)^4
# NCX_scale = C_m/C_myo #Previous descriptors of NCX were whole-cell currents at for cells of ~80pF, adjusting for current in pA
# # 0.077472 : TODO : what is this ? (by Kei)

# """Constants related to the sodium-hydrogen exchanger"""
# n_H = 1
# m_H = 3
# K_H_i_mod = 3.07e-5
# K_H_o_mod = 4.8e-7
# k1_p = 10.5
# k1_m = 0.201
# k2_p = 15.8
# k2_m = 183
# K_Na_i = 16.2
# K_Na_o = 195
# K_H_i = 6.05e-4
# K_H_o = 1.62e-3
# N_NaH_channel = 4899
# I_NaH_scale = 1.00 # for temp 23C
# # I_NaH_scale = 1.29 # for temp 37C

# """Constants related to the calcium pump"""
# I_Ca_ATP_bar = 6.0 #ICaP = 4.0 (pA), Nygren1998
# I_Ca_ATP_scale = 1.0 # for temp 23C
# # I_Ca_ATP_scale = 1.29 # for temp 37C
# k_Ca_ATP = 0.0002 #kCaP = 0.0002 (mmol/L)

# """Constants related to the ultra-rapidly rectifying potassium channel
# g_K_ur = 0.245, gsus = 2.75 nS in Nygren1998"""
# g_K_ur = 1.4623


# """Constants related to I_K_DR"""
# g_K_DR = 0.0289*C_m #constant given in nS/pF per Bob Clark communication 3/13/16, assuming C_m as above
# #g_K_DR = 0.03*C_m #testing
# i_K_DR = 1.0 #per Bob Clark communication 3/13/16, unless we want to remove inactivation 
# act_DR_shift = 10.0 # (mV)

# """Constants related to the two-pore potassium channel"""
# P_K = 3.1e-6*sqrt(5/140) # Not currently used TODO : used ? by Kei 
# #P_K = 10e-6*sqrt(5/140)
# Q = 0.1 #Based on experimental data from Bob; slope of IV in isopotential recording conditions
# I_K_2pore_0 = 0.0
# #I_K_2pore_0 = 1.5
# I_K_2pore_scale = 1.35 # for temp 23C
# # I_K_2pore_scale = 1.29*1.35 # for temp 37C

# """Constants related to the calcium-activated potassium channel"""
# Zj = 0.70
# Vhj = 250
# ZL = 0.1
# L0 = 12e-6
# KDc = 3e-6
# C = 8
# D = 25
# E = 2.4
# Gmax = 3.8*2.4
# N_channel = 1.0
# E_K_Ca_act = 42
# # New Sun, et al formulation
# gBK = 2.50

# """Constants related to the TRP1 channel"""
# g_TRP1 = 1.e-4*0.75*10/4
# a_TRP1 = 80
# b_TRP1 = -1000

# """Constants related to the TRPV4 channel"""
# g_TRPV4 = 0.00046875*0.05 # (nS)
# #g_TRPV4 = 1.0e-1*1.05*100/4
# a_TRPV4 = 80
# b_TRPV4 = -1000
