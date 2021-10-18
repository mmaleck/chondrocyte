"""
This file is part of the chondrocyte modelling project at Simula
Research Laboratory, Norway. Refer to the files README and LICENSE for
more information about the project as well as terms of distribution.
 
Author : Kei Yamamoto, Sofie Fischer, M.M.Malacker
Data created : July, 2021
Python version : 3.8.2
"""

import numpy as np
from math import pi, log, ceil, exp, sqrt
from scipy.signal import square
from scipy.linalg import null_space
from params import params_dict

def rhs(y, t, params_dict):
    V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal = y

    # FIXME: Fixing the volume while debugging 
    vol_i = params_dict["vol_i_0"]

    if (params_dict["apply_Vm"] == True):
        V = appliedVoltage(t)
    else:
        V = V

    #Define external concentrations
    K_o = appliedPotassiumConcentration(t)

    #Calculate background currents
    g_K_b_bar = params_dict["g_K_b_bar"]
    I_Na_b = backgroundSodium(V, Na_i, None, enable_I_Na_b=True)
    I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar, None, enable_I_K_b=True)
    I_Cl_b = backgroundChloride(V, Cl_i, enable_I_Cl_b=True)
    I_leak = backgroundLeak(V, enable_I_leak=False)

    #Calculate pump and exchanger currents
    Na_i_0 = params_dict["Na_i_0"]
    I_NaK = sodiumPotassiumPump(V, K_o, Na_i_0, enable_I_NaK=True)
    I_NaCa = sodiumCalciumExchanger(V, Ca_i, Na_i_0, enable_I_NaCa=True)
    I_NaH = sodiumHydrogenExchanger(Na_i, H_i, enable_I_NaH=True)
    I_Ca_ATP = calciumPump(Ca_i, enable_I_Ca_ATP=True)

    # Calculate potassium currents
    P_K = params_dict["P_K"]; temp = params_dict["temp"]
    I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, enable_I_K_ur=False) # essential reimplemented as I_K_DR, so enable_I_K_ur has to be False
    I_K_DR = DelayedRectifierPotassium(V, enable_I_K_DR=True)
    I_K_2pore = twoPorePotassium(V, K_i, K_o, P_K, enable_I_K_2pore=True)
    I_K_Ca_act = calciumActivatedPotassium(V, Ca_i, enable_I_K_Ca_act=True)
    I_K_ATP = potassiumPump(V, K_i, K_o, None, Na_i, temp, enable_I_K_ATP=True)
  
    # Calculate other currents
    I_ASIC = voltageActivatedHydrogen(enable_I_ASIC=False)
    I_TRP1 = stretchActivatedTrip(V, enable_I_TRP1=False)
    I_TRP2 = osteoArthriticTrip(enable_I_TRP2=False)
    I_TRPV4 = TripCurrent(V, enable_I_TRPV4=False)
    I_stim = externalStimulation(t, enable_I_stim=False)

    #Total ionic contribution (pA)
    I_i = I_Na_b + I_K_b + I_Cl_b + I_leak \
        + I_NaK + I_NaCa + I_Ca_ATP \
        + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP \
        + I_ASIC + I_TRP1 + I_TRP2 + I_TRPV4

    # Determine incremental changes in evolving quantities
    F = params_dict["F"]
    C_m = params_dict["C_m"]

    # Evolve calmodulin - Ca buffer
    cal_dot = 200000.0*Ca_i*(1.0 - cal) - 476.0*cal

    #Evolve the concentrations
    clamp_Na_i = params_dict["clamp_Na_i"]
    clamp_K_i = params_dict["clamp_K_i"]
    calmp_Ca_i = params_dict["calmp_Ca_i"]
    clamp_H_i = params_dict["clamp_H_i"]
    clamp_Cl_i = params_dict["clamp_Cl_i"]

    if (clamp_Na_i == True):
        Na_i_dot = 0
    else:
        Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH + I_Cl_b + 0.5*I_leak)/(vol_i*F)
    
    if (clamp_K_i == True):
        K_i_dot = 0
    else:
        K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP + 0.5*I_leak)/(vol_i*F)
    
    if (calmp_Ca_i == True):
        Ca_i_dot = 0
    else:
        Ca_i_dot =   -(I_Ca_ATP - 2*I_NaCa + I_TRPV4)/(2*vol_i*F) - 0.045*cal_dot
    
    if (clamp_H_i == True):
        H_i_dot = 0
    else:
        H_i_dot =  - (I_NaH)/(vol_i*F)

    if (clamp_Cl_i == True):
        Cl_i_dot = 0
    else:    
        Cl_i_dot =  (I_Cl_b)/(vol_i*F)

    a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur = ultraRapidlyRectifyingPotassiumHelper(V)

    a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur
    i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur

    Na_o = params_dict["Na_o"]
    Ca_o = params_dict["Ca_o"]
    H_o = params_dict["H_o"]
    Cl_o = params_dict["Cl_o"]

    K_i_0 = params_dict["K_i_0"]
    Ca_i_0 = params_dict["Ca_i_0"]
    H_i_0 = params_dict["H_i_0"]
    Cl_i_0 = params_dict["Cl_i_0"]

    #Think this is volume from one of the UK papers...check later 3/16/2016
    osm_i = Na_i + K_i + Ca_i + H_i + Cl_i
    osm_o = Na_o + K_o + Ca_o + H_o + Cl_o
    P_f = 10.2e-4                                      # water permeability of the cell (cm sec^-1), check unit 
    SA = 6.0**(2.0/3.0)*pi**(1.0/3.0)*vol_i**(2.0/3.0) # Surface Area
    V_W = 18.0                                         # Molar volume of water, check unit  
    # compute dvolume/dt 
    vol_i_dot = P_f*SA*V_W*(osm_i - osm_o)*1e-4

    # print(vol_i_dot)

    apply_Vm = params_dict["apply_Vm"]
    if (apply_Vm == True):
        V_dot = 0.0
    else:
        V_dot = 1/C_m*(-I_i + I_stim)

    clamp_conc = params_dict["clamp_conc"]
    if (clamp_conc == True):
        Na_i_dot = 0.0
        K_i_dot = 0.0
        Ca_i_dot= 0.0
        H_i_dot = 0.0
        Cl_i_dot = 0.0
    
    # Return RHS as sequence
    return (V_dot, Na_i_dot, K_i_dot, Ca_i_dot, H_i_dot, Cl_i_dot, a_ur_dot, i_ur_dot, vol_i_dot, cal_dot)

def nernstPotential(z, X_i, X_o):
    """Potential of an ion X across the membrane (mV)."""
    R = params_dict["R"]; T = params_dict["T"]; F = params_dict["F"]
    E_X = (R*T)/(z*F)*log(X_o/X_i)
    return E_X

def appliedVoltage(t):
    """Applied voltage (mV)"""
    clamp_Vm = params_dict["clamp_Vm"]; ramp_Vm = params_dict["ramp_Vm"]; step_Vm = params_dict["step_Vm"]
    if (clamp_Vm == True):
        V_0 = params_dict["V_0"]
        V = V_0
    elif (ramp_Vm == True):
        V_0 = params_dict["V_0"]; V_final = params_dict["V_final"]; t_final = params_dict["t_final"]
        V = V_0 + (V_final - V_0)*t/t_final
    elif (step_Vm == True):
        t_cycle = params_dict["t_cycle"]; t_stim = params_dict["t_stim"]
        V = (ceil((t - 30)/t_cycle)*square((t - 30)*2*pi/t_cycle, t_stim/t_cycle) + ceil((t - 30)/t_cycle))/2*10 - 90

    if (V == 0):
        V = 0.01

    return V

def appliedPotassiumConcentration(t):
    step_K_o = params_dict["step_K_o"]; K_o_0 = params_dict["K_o_0"]
    if (step_K_o == False):
        K_o = K_o_0
    else:
        if (t <= 10):
            K_o = 5
        elif (t > 10 & t <= 20):
            K_o = 30
        elif (t > 20 & t <= 30):
            K_o = 75
        elif (t > 30 & t <= 40):
            K_o = 140
        else:
            K_o = 5

    return K_o

def backgroundSodium(V, Na_i, E_Na, enable_I_Na_b):
    """Background sodium current from "Ionic channels of excitable
    membranes," B. Hille. (pA)
    """
    if (enable_I_Na_b == True):
        z_Na = params_dict["z_Na"]; g_Na_b_bar = params_dict["g_Na_b_bar"]; Na_o = params_dict["Na_o"]
        I_Na_b_scale = params_dict["I_Na_b_scale"]
        if E_Na == None:
            E_Na = nernstPotential(z_Na, Na_i, Na_o)
        I_Na_b = I_Na_b_scale*g_Na_b_bar*(V - E_Na)
    else:
        I_Na_b = 0.0

    return I_Na_b

def backgroundPotassium(V, K_i, K_o, g_K_b_bar, E_K, enable_I_K_b):
    """Background potassium current from "Ionic channels of excitable
    membranes," B. Hille. (pA)
    """
    if (enable_I_K_b == True):
        z_K = params_dict["z_K"]
        if E_K == None :
            E_K = nernstPotential(z_K, K_i, K_o)
        I_K_b = g_K_b_bar*(V - E_K)
    else:
        I_K_b = 0.0

    return I_K_b

def backgroundChloride(V, Cl_i, enable_I_Cl_b):
    """Background chloride current from "Ionic channels of excitable
    membranes," B. Hille. (pA)
    """
    if (enable_I_Cl_b == True):
        z_Cl = params_dict["z_Cl"]; g_Cl_b_bar = params_dict["g_Cl_b_bar"]; Cl_o = params_dict["Cl_o"]
        #E_Cl = nernstPotential(z_Cl, Cl_o, Cl_i)
        #E_Cl = -40.0
        E_Cl = params_dict["E_Cl"]
        I_Cl_b = g_Cl_b_bar*(V - E_Cl)
    else:
        I_Cl_b = 0.0

    return I_Cl_b

def backgroundLeak(V, enable_I_leak):
    if (enable_I_leak == True):
        g_leak = params_dict["g_leak"]
        I_leak = g_leak*V
    else:
        I_leak = 0.0
    
    return I_leak


def sodiumPotassiumPump(V, K_o, Na_i_0, enable_I_NaK):
    """Sodium-potassium pump from "Mathematical Model of an Adult Human
    Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
    Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
    Giles. Circ. Res. 1998; 82; 63-81 (Table 12, pp. 77) (pA)
    """
    if (enable_I_NaK == True):
        I_NaK_bar = params_dict["I_NaK_bar"]; K_NaK_K = params_dict["K_NaK_K"]; K_NaK_Na = params_dict["K_NaK_Na"]
        I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) \
                *(Na_i_0**(1.5)/(Na_i_0**(1.5) + K_NaK_Na**(1.5))) \
                *(V + 150.0)/(V + 200.0)
    else:
        I_NaK = 0.0

    return I_NaK

def sodiumCalciumExchanger(V, Ca_i, Na_i_0, enable_I_NaCa):
    """Sodium-calcium exchanger from "Mathematical Model of an Adult Human
    Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
    Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
    Giles. Circ. Res. 1998; 82; 63-81 (Table 13, pp. 77) (pA)
    """
    if (enable_I_NaCa == True):
        F = params_dict["F"]; R = params_dict["R"]; T = params_dict["T"]
        Na_o = params_dict["Na_o"]; Ca_o = params_dict["Ca_o"]
        K_NaCa = params_dict["K_NaCa"]; gamma_Na = params_dict["gamma_Na"]; d_NaCa = params_dict["d_NaCa"]
        NCX_scale = params_dict["NCX_scale"]

        I_NaCa = NCX_scale*K_NaCa*(Na_i_0**3*Ca_o*exp(gamma_Na*V*F/(R*T)) \
                - Na_o**3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
                /(1.0 + d_NaCa*(Na_o**3*Ca_i + Na_i_0**3*Ca_o)) 
    else:
        I_NaCa = 0.0
    
    return I_NaCa

def sodiumHydrogenExchanger(Na_i, H_i, enable_I_NaH):
    """Sodium-hydrogen exchanger from "A Model of Na+/H+ Exchanger and Its
    Central Role in Regulation of pH and Na+ in Cardiac Myocytes," Chae
    Young Cha, Chiaki Oka, Yung E. Earm, Shigeo Wakabayashi, and Akinori
    Noma. Biophysical Journal 2009; 97; 2674-2683 (pp. 2675) (pA)
    """
    if (enable_I_NaH == True):
        n_H = params_dict["n_H"]; K_H_i_mod = params_dict["K_H_i_mod"]; I_NaH_scale = params_dict["I_NaH_scale"]
        k1_p = params_dict["k1_p"]; k1_m = params_dict["k1_m"]; k2_p = params_dict["k2_p"]; k2_m = params_dict["k2_m"]
        Na_o = params_dict["Na_o"]; H_o = params_dict["H_o"]; N_NaH_channel = params_dict["N_NaH_channel"]
        K_Na_o = params_dict["K_Na_o"]; K_H_o = params_dict["K_H_o"];K_Na_i = params_dict["K_Na_i"]; K_H_i = params_dict["K_H_i"]

        I_NaH_mod  = 1/(1 + (K_H_i_mod**(n_H)/H_i**(n_H)))
        t1 = k1_p*Na_o/K_Na_o / (1 + Na_o/K_Na_o + H_o/K_H_o)
        t2 = k2_p*H_i/K_H_i   / (1 + Na_i/K_Na_i + H_i/K_H_i)
        t3 = k1_m*Na_i/K_Na_i / (1 + Na_i/K_Na_i + H_i/K_H_i)
        t4 = k2_m*H_o/K_H_o   / (1 + Na_o/K_Na_o + H_o/K_H_o)
        I_NaH_exch = (t1*t2 - t3*t4) / (t1 + t2 + t3 + t4)
        I_NaH = I_NaH_scale*N_NaH_channel*I_NaH_mod*I_NaH_exch
    else:
        I_NaH = 0.0
    
    return I_NaH


def calciumPump(Ca_i, enable_I_Ca_ATP):
    """Calcium pump from Nygren et al. (pA)"""
    if (enable_I_Ca_ATP == True):
        I_Ca_ATP_bar = params_dict["I_Ca_ATP_bar"]; k_Ca_ATP = params_dict["k_Ca_ATP"]; I_Ca_ATP_scale = params_dict["I_Ca_ATP_scale"]
        I_Ca_ATP = I_Ca_ATP_scale*I_Ca_ATP_bar*(Ca_i/(Ca_i + k_Ca_ATP))
    else:
        I_Ca_ATP = 0.0

    return I_Ca_ATP

def ultraRapidlyRectifyingPotassiumHelper(V):
    """Ultra-rapidly rectifying potassium channel from "Action potential rate
    dependence in the human atrial myocyte," M. M. Maleckar, J. L.
    Greenstein, W. R. Giles and N. A. Trayanova. Am. J. Physiol. Heart.
    Circ. Physiol. 2009; 297; 1398-1410 (Appendix, pp. 1408) (pA)
    """
    a_ur_inf   = 1.0/(1.0 + exp(-(V + 26.7)/4.1))
    i_ur_inf   = 1.0/(1.0 + exp((V - 30.0)/10.0))
    tau_a_ur   = 0.005/(1.0 + exp((V + 5.0)/12.0))
    tau_i_ur   = 0.59/(1.0 + exp((V + 10.0)/24.0)) + 0.01

    return a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur

def ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, enable_I_K_ur):
    if (enable_I_K_ur == True):
        z_K = params_dict["z_K"]; g_K_ur = params_dict["g_K_ur"]
        a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur = ultraRapidlyRectifyingPotassiumHelper(V)
        E_K        = nernstPotential(z_K, K_i, K_o)
        I_K_ur     = g_K_ur*a_ur*(V - E_K)
    else:
        I_K_ur = 0.0
    return I_K_ur
   
def DelayedRectifierPotassium(V, enable_I_K_DR):
    """ Delayed rectifer potassium channel from Clark, et al (pA)"""
    if (enable_I_K_DR == True):
        g_K_DR = params_dict["g_K_DR"]; i_K_DR = params_dict["i_K_DR"]; act_DR_shift = params_dict["act_DR_shift"]
        alpha_K_DR = 1.0/(1.0 + exp(-(V + 26.7 + act_DR_shift)/4.1)) # Same as what I had. 3/3/16
        # E_K        = nernstPotential(z_K, K_i, K_o)
        E_K = -83
        I_K_DR = g_K_DR*alpha_K_DR*i_K_DR*(V - E_K)
    else:
        I_K_DR = 0.0
    
    return I_K_DR

def ultrarapidlyRectifyingPotassium_ref(V, K_i, K_o):
    """ From Bob Clark et al., J. Physiol. 2011, Figure 4 - IKDR"""
    G_K = params_dict["G_K"]
    V_h = params_dict["V_h"]
    S_h = params_dict["S_h"]
    C_m = params_dict["C_m"]; z_K = params_dict["z_K"]
    E_K        = nernstPotential(z_K, K_i, K_o)
    I_K_ur_ref = G_K*(V - E_K)/(1 + exp(-(V - V_h)/S_h)) * C_m/1000.0

    return I_K_ur_ref

def twoPorePotassium(V, K_i_0, K_o, Q, enable_I_K_2pore):
    """ Two-pore potassium current - novel?; modeled as a simple Boltzmann
    relationship via GHK (pA)"""
    if (enable_I_K_2pore == True):
        F = params_dict["F"]; R = params_dict["R"]; T = params_dict["T"]; z_K = params_dict["z_K"]; I_K_2pore_0 = params_dict["I_K_2pore_0"]
        I_K_2pore_scale = params_dict["I_K_2pore_scale"]
        if V == 0: # based on Flax, Matt R. and Holmes, W.Harvey (2008) Goldman-Hodgkin-Katz Cochlear Hair Cell Models - a Foundation for Nonlinear Cochlear Mechanics, Conference proceedings: Interspeech 2008
            I_K_2pore = I_K_2pore_scale*5*Q*sqrt(K_o/K_i_0)*R*T/(z_K*F)*(1-(K_o/K_i_0)) + I_K_2pore_0
        else:
            I_K_2pore = I_K_2pore_scale*5*Q*sqrt(K_o/K_i_0)*V*(1 - (K_o/K_i_0)*exp(-z_K*V*F/(R*T)))/(1- exp(-z_K*V*F/(R*T))) + I_K_2pore_0
    else:
        I_K_2pore = 0.0

    return I_K_2pore

# FIXME: Check the following carefully
def calciumActivatedPotassium(V, Ca_i, enable_I_K_Ca_act):
    """Calcium-activated potassium current - Sun, et al formulations (pA)"""
    if (enable_I_K_Ca_act == True):
        F = params_dict["F"]; R = params_dict["R"]; T = params_dict["T"]
        # I_K_Ca_act (new version) (pA), with converted Ca_i units for model
        # Set constants
        convert_units = 1e6 # Convert from nM (e-9) to mM (e-3)
        gBK = params_dict["gBK"]
        E_K = -83  # Sun, et al
        #tspan = 1e-3*[0:1:300] %time in seconds
        #C_m = 7 %pF

        K_C = 17
        K_O = 0.5

        A = np.array([0.659,3.955,25.05,129.2,261.1])*(1/65)
        B = np.array([2651.7, 1767.8, 1244.0, 713.0, 160.0])*(4.5)
        
        z_CO = 0.718
        z_OC = 0.646

        FVonRT = F*V/(R*T)
        
        alpha_BK = np.zeros(5)
        beta_BK = np.zeros(5)

        for m in range(5):
            alpha_BK[m] = A[m]*exp(+z_CO*FVonRT)
            beta_BK[m]  = B[m]*exp(-z_OC*FVonRT)
            #delta(i) = C(i)*exp(+z_OI*FVonRT)
            #gamma(i) = D(i)*exp(-z_IO*FVonRT)

        # Build Markov Matrix
        n = 10
        M = np.zeros((n,n))

        # numbering scheme, offsets:  
        closed = 0
        open = 5
        #inactivated = 10

        # Vertical transitions
        for k in range(5):
            M[closed + k , open + k] = alpha_BK[k]
            M[open + k , closed + k] = beta_BK[k]

            #M(open + i, inactivated + i) = delta(i)
            #M(inactivated + i, open + i) = gamma(i)
    
        #  Horizontal transitions
        for jj in range(4):
            # on rates:
            k_on = (4-jj)*Ca_i*convert_units

            M[closed + jj, closed + jj + 1] = k_on
            M[open + jj, open + jj + 1] = k_on
           #M(inactivated + i, inactivated + i + 1) = k_on

        # off rates:
            M[closed + jj + 1, closed + jj] = (jj+1)*K_C
            M[open + jj + 1, open + jj] = (jj+1)*K_O
        #  M(inactivated + i + 1, inactivated + i) = i*K_I
        
    # Transpose since we have used above: M(from, to) = k_{from, to}
        M = np.transpose(M)
        for kk in range(n):
            M[kk,kk] = -np.sum(M[:,kk], axis=0)

    #  Solve the system for Steady state at Ca_i_ss and V     
        eq = null_space(M) #find nullspace for BK (equilibrium)
        eq = eq/np.sum(eq) #Find unique equilibrium by scaling probabilities
        #ode = @(t,y) ode_system(t,y,V)
        #[T,S] = ode15s(ode,tspan,eq)
        open = np.sum(eq[5:10]) #Calculate total steady-state open probability
        I_BK = gBK*open*(V-E_K) #Calculate steady-state current in pA
        I_K_Ca_act = I_BK
        
    else:
        I_K_Ca_act = 0.0
    
    return I_K_Ca_act
    
def potassiumPump(V, K_i, K_o, E_K, Na_i, temp, enable_I_K_ATP):
    """ATP-dependent K+ current, based on Simulation of Action Potentials From Metabolically Impaired Cardiac Myocytes Role of ATP-Sensitive K+ Current by Ferreo et al, 1996"""
    if (enable_I_K_ATP == True):
        F = params_dict["F"]; R = params_dict["R"]; 
        T = temp + 273.15
        # compute f_N
        K_h_Na_0 = params_dict["K_h_Na_0"]
        delta_Na = params_dict["delta_Na"]
        gamma_zero = 33.375*(K_o/5.4)**(0.24)/500 
        K_h_Na = K_h_Na_0*exp(-(delta_Na*F*V)/(R*T))
        f_N = 1/(1+(Na_i/K_h_Na)**2)
        # compute f_M
        delta_Mg = params_dict["delta_Na"]; Mg_i = params_dict["Mg_i"]
        K_h_Mg = 0.65/sqrt(K_o+5)*exp(-2*delta_Mg*F*V/(R*T))
        f_M = 1/(1+(Mg_i/K_h_Mg))
        # compute f_T
        Q_10 = params_dict["Q_10"]
        f_T = Q_10**((temp-36)/10)
        # compute g_0 based on f_M, f_N, and f_T
        g_0 = gamma_zero*f_M*f_N*f_T
        # compute f_ATP
        H_K_ATP = params_dict["H_K_ATP"]
        K_m_ATP = params_dict["K_m_ATP"]
        C_A = params_dict["C_A"]
        V_0 = params_dict["V_0"]
        ATP_i = params_dict["ATP_i"]
        ADP_i = C_A - ATP_i
        H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i)
        K_m = 35.8 + 17.9*ADP_i**(K_m_ATP)
        f_ATP = 1.0/(1.0 + (ATP_i/K_m)**H)

        z_K = params_dict["z_K"]
        # if E_K == None:
            # E_K = nernstPotential(z_K, K_i, K_o)
        # from IPython import embed; embed(); exit(1)
        E_K = nernstPotential(z_K, K_i, K_o)
        sigma = params_dict["sigma"]; p_0 = params_dict["p_0"] 
        I_K_ATP = sigma*g_0*p_0*f_ATP*(V - E_K)
    else:
        I_K_ATP = 0.0
    
    return I_K_ATP

def externalStimulation(t, enable_I_stim):
    """External stimulation"""
    if (enable_I_stim == True):
        t_cycle = params_dict["t_cycle"]; t_stim = params_dict["t_stim"]; I_stim_bar = params_dict["I_stim_bar"]
        I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle)
    else :
        I_stim = 0.0
    
    return I_stim

# FIXME: Implement the voltage-activated hydrogen channel
def voltageActivatedHydrogen(enable_I_ASIC):
    if (enable_I_ASIC == True):
        I_ASIC = 0.0
    else:
        I_ASIC = 0.0
    
    return I_ASIC

def TripCurrent(V, enable_I_TRPV4):
    """Implement the TRPV4 channel"""
    if (enable_I_TRPV4 == True):
        g_TRPV4 = params_dict["g_TRPV4"];  a_TRPV4 = params_dict["a_TRPV4"]; b_TRPV4 = params_dict["b_TRPV4"]
        if(V < 0):
            I_TRPV4 = g_TRPV4*(b_TRPV4*V + (1 - b_TRPV4)*a_TRPV4*(1 - (1 - (V/a_TRPV4))*(1 - (V/a_TRPV4))*(1-(V/a_TRPV4))))
        else:
            I_TRPV4 = 2*g_TRPV4*V**3
    else:
        I_TRPV4 = 0.0

    return I_TRPV4

# FIXME: Implement the stretch-activated TRP channel    
def stretchActivatedTrip(V, enable_I_TRP1):
    if (enable_I_TRP1 == True):
        g_TRP1 = params_dict["g_TRP1"]; a_TRP1 = params_dict["a_TRP1"]; b_TRP1 = params_dict["b_TRP1"]
        if(V < 0):
            I_TRP1 = g_TRP1*(b_TRP1*V + (1 - b_TRP1)*a_TRP1*(1 - (1 - (V/a_TRP1))*(1 - (V/a_TRP1))*(1-(V/a_TRP1))))
        else:
            I_TRP1 = 2*g_TRP1*V**3
    else:
        I_TRP1 = 0.0

    return I_TRP1

# FIXME: Implement the osteo-arthritic TRP channel 
def osteoArthriticTrip(enable_I_TRP2):
    if (enable_I_TRP2 == True):
        I_TRP2 = 0.0
    else:
        I_TRP2 = 0.0

    return I_TRP2