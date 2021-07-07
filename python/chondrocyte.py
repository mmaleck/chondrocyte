import params
import numpy as np
from math import pi, log, ceil, exp, sqrt
from scipy.signal import square
from scipy.linalg import null_space

def rhs(t, y, g_K_b_bar, P_K, Gmax):
    V, Na_i, K_i, Ca_i, H_i, Cl_i, a_ur, i_ur, vol_i, cal = y

    if (params.apply_Vm == True):
        V = appliedVoltage(t)
    else:
        V = V

    
    #Define external concentrations
    K_o = appliedPotassiumConcentration(t)

    #Calculate background currents
    I_Na_b = backgroundSodium(V, Na_i)
    I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar)
    I_Cl_b = backgroundChloride(V, Cl_i)
    I_leak = backgroundLeak(V)

     #Calculate pump and exchanger currents
    I_NaK = sodiumPotassiumPump(V, Na_i, K_i, K_o)
    I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i)
    I_NaH = sodiumHydrogenExchanger(Na_i, H_i)
    I_Ca_ATP = calciumPump(Ca_i)

    # Calculate potassium currents
    I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur)
    I_K_DR = DelayedRectifierPotassium(V)
    I_K_2pore = twoPorePotassium(V, K_i, K_o, P_K)
    I_K_Ca_act = calciumActivatedPotassium(V, Ca_i)
    I_K_ATP = potassiumPump(V, K_i, K_o)
  
    # Calculate other currents
    I_ASIC = voltageActivatedHydrogen()
    I_TRP1 = stretchActivatedTrip(V)
    I_TRP2 = osteoArthriticTrip()
    I_TRPV4 = TripCurrent(V)
    I_stim = externalStimulation(t)

    #Total ionic contribution (pA)
    I_i = I_Na_b + I_K_b + I_Cl_b + I_leak \
      + I_NaK + I_NaCa + I_Ca_ATP \
      + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP \
      + I_ASIC + I_TRP1 + I_TRP2 + I_TRPV4

    # Determine incremental changes in evolving quantities
    F = params.F
    C_m = params.C_m

    # Evolve calmodulin - Ca buffer
    cal_dot = 200000.0*Ca_i*(1.0 - cal) - 476.0*cal
    # cal_dot = 200000.0*Ca_i*(1.0 - cal);

    #Evolve the concentrations
    clamp_Na_i = params.clamp_Na_i
    clamp_K_i = params.clamp_K_i
    if (clamp_Na_i == True):
        Na_i_dot = 0
    else:
        Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH + I_Cl_b + 0.5*I_leak)/(vol_i*F)
        #Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH)/(vol_i*F)
    
    if (clamp_K_i == True):
        K_i_dot = 0
    else:
        K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP + 0.5*I_leak)/(vol_i*F)
        #K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP)/(vol_i*F)
    
    # Ca_i_dot =   -(I_Ca_ATP - 2*I_NaCa + I_TRPV4)/(2*vol_i*F) - 0.045*cal_dot
    Ca_i_dot =   -(I_Ca_ATP - 2*I_NaCa + I_TRPV4)/(2*vol_i*F) - 0.045*cal_dot
    #H_i_dot =  - (I_NaH)/(vol_i*F)
    H_i_dot = 0
    Cl_i_dot =  (I_Cl_b)/(vol_i*F)
    #Cl_i_dot = 0

    a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur = ultraRapidlyRectifyingPotassiumHelper(V)

    a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur
    i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur


    Na_o = params.Na_o
    Ca_o = params.Ca_o
    H_o = params.H_o
    Cl_o = params.Cl_o

    Na_i_0 = params.Na_i_0
    K_i_0 = params.K_i_0
    Ca_i_0 = params.Ca_i_0
    H_i_0 = params.H_i_0
    Cl_i_0 = params.Cl_i_0
    Na_i_clamp = params.Na_i_clamp


   #Think this is volume from one of the UK papers...check later 3/16/2016
    osm_i_0 = Na_i_0 + K_i_0 + Ca_i_0 + H_i_0 + Cl_i_0
    osm_i = Na_i + K_i + Ca_i + H_i + Cl_i
    osm_o = Na_o + K_o + Ca_o + H_o + Cl_o
    dosm = osm_i_0 - osm_o

    P_f = 10.2e-4
    SA = 6.0^(2.0/3.0)*pi^(1.0/3.0)*vol_i^(2.0/3.0)
    V_W = 18.0
    vol_i_dot = P_f*SA*V_W*(osm_i - osm_o - dosm)

    apply_Vm = params.apply_Vm
    if (apply_Vm == True):
        V_dot = 0.0
    else:
        V_dot = 1/C_m*(-I_i + I_stim)

    clamp_conc = params.clamp_conc
    if (clamp_conc == True):
        Na_i_dot = 0.0
        K_i_dot = 0.0
        Ca_i_dot= 0.0
        H_i_dot = 0.0
        Cl_i_dot = 0.0
    
    # Return RHS as sequence
    return (V_dot, Na_i_dot, K_i_dot, Ca_i_dot, H_i_dot, Cl_i_dot, a_ur_dot, i_ur_dot, vol_i_dot, cal_dot)


# Potential of an ion X across the membrane (mV).
def nernstPotential(z, X_i, X_o):
  R = params.R, T = params.T, F = params.F
  E_X = (R*T)/(z*F)*log(X_o/X_i)
  return E_X

# Applied voltage (mV)
def appliedVoltage(t):
    clamp_Vm = params.clamp_Vm, ramp_Vm = params.ramp_Vm, step_Vm = params.step_Vm
    if (clamp_Vm == True):
        V_0 = params.V_0
        V = V_0
    elif (ramp_Vm == True):
        V_0 = params.V_0, V_final = params.V_final, t_final = params.t_final
        V = V_0 + (V_final - V_0)*t/t_final
    elif (step_Vm == True):
        t_cycle = params.t_cycle, t_stim = params.t_stim
        # the following needs to be verified 
        V = (ceil((t - 30)/t_cycle)*square((t - 30)*2*pi/t_cycle, t_stim/t_cycle) + ceil((t - 30)/t_cycle))/2*10 - 90

    if (V == 0):
        V = 0.01

    return V

def appliedPotassiumConcentration(t):
    step_K_o = params.step_K_O, K_o_0 = params.K_o_0
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

# Background sodium current from "Ionic channels of excitable
# membranes," B. Hille. (pA)
def backgroundSodium(V, Na_i):
    enable_I_Na_b = params.enable_I_Na_b
    if (enable_I_Na_b == True):
        z_Na = params.z_Na, g_Na_b_bar = params.g_Na_b_bar, Na_o = params.Na_o
        E_Na = nernstPotential(z_Na, Na_i, Na_o)
        I_Na_b = g_Na_b_bar*(V - E_Na)
    else:
        I_Na_b = 0.0

    return I_Na_b

# Background potassium current from "Ionic channels of excitable
# membranes," B. Hille. (pA)
def backgroundPotassium(V, K_i, K_o, g_K_b_bar):
    enable_I_K_b = params.enable_I_K_b
    if (enable_I_K_b == True):
        z_K = params.z_K
        E_K = nernstPotential(z_K, K_i, K_o)
        I_K_b = g_K_b_bar*(V - E_K)
    else:
        I_K_b = 0.0

    return I_K_b

# Background chloride current from "Ionic channels of excitable
# membranes," B. Hille. (pA)
def backgroundChloride(V, Cl_i):
    enable_I_Cl_b = params.enable_I_Cl_b
    if (enable_I_Cl_b == True):
        z_Cl = params.z_Cl, g_Cl_b_bar = params.g_Cl_b_bar, Cl_o = params.Cl_o
        #E_Cl = nernstPotential(z_Cl, Cl_o, Cl_i)
        #E_Cl = -40.0
        E_Cl = -65.0
        I_Cl_b = g_Cl_b_bar*(V - E_Cl)
    else:
        I_Cl_b = 0.0

    return I_Cl_b

def backgroundLeak(V):
    enable_I_leak = params.enable_I_leak
    if (enable_I_leak == True):
        g_leak = params.g_leak
        I_leak = g_leak*V
    else:
        I_leak = 0.0
    
    return I_leak

'''
Sodium-potassium pump from "Mathematical Model of an Adult Human
Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
Giles. Circ. Res. 1998; 82; 63-81 (Table 12, pp. 77) (pA)
'''
def sodiumPotassiumPump(V, Na_i, K_i, K_o):
    enable_I_NaK = params.enable_I_NaK , Na_i_0 = params.Na_i_0
    if (enable_I_NaK == True):
        I_NaK_bar = params.I_NaK_bar, K_NaK_K = params.K_NaK_K, K_NaK_Na = params.K_NaK_Na
        I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) \
            *(Na_i_0**(1.5)/(Na_i_0**(1.5) + K_NaK_Na**(1.5))) \
            *(V + 150.0)/(V + 200.0)
    else:
        I_NaK = 0.0

    return I_NaK

'''
Sodium-calcium exchanger from "Mathematical Model of an Adult Human
Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
Giles. Circ. Res. 1998; 82; 63-81 (Table 13, pp. 77) (pA)
'''
def sodiumCalciumExchanger(V, Na_i, Ca_i):
    enable_I_NaCa = params.enable_I_NaCa
    if (enable_I_NaCa == True):
        F = params.F, R = params.R, T = params.T
        Na_o = params.Na_o, Ca_o = params.Ca_o
        K_NaCa = params.K_NaCa, gamma_Na = params.gamma_Na, d_NaCa = params.d_NaCa, Na_i_0 = params.Na_i_0
        NCX_scale = params.NCX_scale

        I_NaCa = NCX_scale*K_NaCa*(Na_i_0**3*Ca_o*exp(gamma_Na*V*F/(R*T)) \
                        - Na_o**3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) \
                /(1.0 + d_NaCa*(Na_o**3*Ca_i + Na_i_0**3*Ca_o)) 
    else:
        I_NaCa = 0.0
    
    return I_NaCa

'''
Sodium-hydrogen exchanger from "A Model of Na+/H+ Exchanger and Its
Central Role in Regulation of pH and Na+ in Cardiac Myocytes," Chae
Young Cha, Chiaki Oka, Yung E. Earm, Shigeo Wakabayashi, and Akinori
Noma. Biophysical Journal 2009; 97; 2674-2683 (pp. 2675) (pA)
'''
def sodiumHydrogenExchanger(Na_i, H_i):
    enable_I_NaH = params.enable_I_NaH
    if (enable_I_NaH == True):
        n_H = params.n_H, K_H_i_mod = params.K_H_i_mod
        k1_p = params.k1_p, k1_m = params.k1_m, k2_p = params.k2_p, k2_m = params.k2_m
        Na_o = params.Na_o, H_o = params.H_o, N_NaH_channel = params.N_NaH_channel
        K_Na_o = params.K_Na_o, K_H_o = params.K_H_o, K_Na_i = params.K_Na_i, K_H_i = params.K_H_i

        I_NaH_mod  = 1/(1 + (K_H_i_mod^n_H/H_i^n_H))
        t1 = k1_p*Na_o/K_Na_o / (1 + Na_o/K_Na_o + H_o/K_H_o)
        t2 = k2_p*H_i/K_H_i   / (1 + Na_i/K_Na_i + H_i/K_H_i)
        t3 = k1_m*Na_i/K_Na_i / (1 + Na_i/K_Na_i + H_i/K_H_i)
        t4 = k2_m*H_o/K_H_o   / (1 + Na_o/K_Na_o + H_o/K_H_o)
        I_NaH_exch = (t1*t2 - t3*t4) / (t1 + t2 + t3 + t4)
        I_NaH = N_NaH_channel*I_NaH_mod*I_NaH_exch
    else:
        I_NaH = 0.0
    
    return I_NaH

# Calcium pump from Nygren et al. (pA)
def calciumPump(Ca_i):
    enable_I_Ca_ATP = params.enable_I_Ca_ATP
    if (enable_I_Ca_ATP == True):
        I_Ca_ATP_bar = params.I_Ca_ATP_bar, k_Ca_ATP = params.k_Ca_ATP
        I_Ca_ATP = I_Ca_ATP_bar*(Ca_i/(Ca_i + k_Ca_ATP))
    else:
        I_Ca_ATP = 0.0

    return I_Ca_ATP

'''
Ultra-rapidly rectifying potassium channel from "Action potential rate
dependence in the human atrial myocyte," M. M. Maleckar, J. L.
Greenstein, W. R. Giles and N. A. Trayanova. Am. J. Physiol. Heart.
Circ. Physiol. 2009; 297; 1398-1410 (Appendix, pp. 1408) (pA)
'''

def ultraRapidlyRectifyingPotassiumHelper(V):
    a_ur_inf   = 1.0/(1.0 + exp(-(V + 26.7)/4.1))
    i_ur_inf   = 1.0/(1.0 + exp((V - 30.0)/10.0))
    tau_a_ur   = 0.005/(1.0 + exp((V + 5.0)/12.0))
    tau_i_ur   = 0.59/(1.0 + exp((V + 10.0)/24.0)) + 0.01

    return a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur

def ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur):
    enable_I_K_ur = params.enable_I_K_ur
    if (enable_I_K_ur == True):
        z_K = params.z_K, g_K_ur = params.g_K_ur
        # why this function is called ?
        a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur = ultraRapidlyRectifyingPotassiumHelper(V)
        E_K        = nernstPotential(z_K, K_i, K_o)
        I_K_ur     = g_K_ur*a_ur*(V - E_K)
    else:
        I_K_ur = 0.0
    return I_K_ur
   
# Delayed rectifer potassium channel from Clark, et al (pA)
def DelayedRectifierPotassium(V):
    enable_I_K_DR = params.enable_I_K_DR
    if (enable_I_K_DR == True):
        g_K_DR = params.g_K_DR ,i_K_DR = params.i_K_DR, act_DR_shift = params.act_DR_shift 
        alpha_K_DR = 1.0/(1.0 + exp(-(V + 26.7 + act_DR_shift)/4.1)); # Same as what I had. 3/3/16
        # E_K        = nernstPotential(z_K, K_i, K_o)
        E_K = -83
        I_K_DR = g_K_DR*alpha_K_DR*i_K_DR*(V - E_K)
    else:
        I_K_DR = 0.0
    
    return I_K_DR

# From Bob Clark et al., J. Physiol. 2011, Figure 4 - IKDR
def ultrarapidlyRectifyingPotassium_ref(V, K_i, K_o):
    G_K = 28.9;  # pS/pF
    V_h = -26.7; # mV
    S_h = 4.1;   # mV
    C_m = params.C_m, z_K = params.z_K
    E_K        = nernstPotential(z_K, K_i, K_o)
    I_K_ur_ref = G_K*(V - E_K)/(1 + exp(-(V - V_h)/S_h)) * C_m/1000.0

    return I_K_ur_ref

# Two-pore potassium current - novel?; modeled as a simple Boltzmann
# relationship via GHK (pA)
def twoPorePotassium(V, K_i_0, K_o, Q):
    enable_I_K_2pore = params.enable_I_K_2pore
    if (enable_I_K_2pore == True):
        F = params.F, R = params.R, T = params.T, z_K = params.z_K, I_K_2pore_0 = params.I_K_2pore_0
        # OLD, via Harish: I_K_2pore = P_K*z_K^2*V*F^2/(R*T)*(K_i_0 - K_o*exp(-z_K*V*F/(R*T)))/(1- exp(-z_K*V*F/(R*T))) + I_K_2pore_0;
        I_K_2pore = 5*Q*sqrt(K_o/K_i_0)*V*(1 - (K_o/K_i_0)*exp(-z_K*V*F/(R*T)))/(1- exp(-z_K*V*F/(R*T))) + I_K_2pore_0
    else:
        I_K_2pore = 0.0

    return I_K_2pore

'''
Calcium-activated potassium current - Sun, et al formulations (pA)
FIXME: Check the following carefully
'''
def calciumActivatedPotassium(V, Ca_i):
    enable_I_K_Ca_act = params.enable_I_K_Ca_act
    if (enable_I_K_Ca_act == True):
        F = params.F, R = params.R, T = params.T
        # I_K_Ca_act (new version) (pA), with converted Ca_i units for model
        # Set constants
        convert_units = 1e-6 # Convert from nM (e-9) to mM (e-3)
        gBK = 2.50 # Gmax, nS
        E_K = -83  # Sun, et al
        #tspan = 1e-3*[0:1:300] %time in seconds
        #C_m = 7 %pF

        K_C = 17
        K_O = 0.5

        A = 1/65*[0.659,3.955,25.05,129.2,261.1]
        B = 4.5*[2651.7, 1767.8, 1244.0, 713.0, 160.0]

        z_CO = 0.718
        z_OC = 0.646

        FVonRT = F*V/(R*T)
        
        alpha_BK = np.zeros(5)
        beta_BK = np.zeros(5)

        for m in range(5):
            alpha_BK(m) = A(m)*exp(+z_CO*FVonRT)
            beta_BK(m)  = B(m)*exp(-z_OC*FVonRT)
            #delta(i) = C(i)*exp(+z_OI*FVonRT);
            #gamma(i) = D(i)*exp(-z_IO*FVonRT);

        # Build Markov Matrix
        n = 10; 
        M = np.zeros(n,n)

        # numbering scheme, offsets:  
        closed = 0
        open = 5
        #inactivated = 10

        
        # Vertical transitions
        for k in range(5):
            M(closed + k , open + k) = alpha_BK(k)
            M(open + k , closed + k) = beta_BK(k)

            #M(open + i, inactivated + i) = delta(i);
            #M(inactivated + i, open + i) = gamma(i);
    

        #  Horizontal transitions
        for jj in range(4):
            # on rates:
            k_on = (5-jj)*Ca_i*convert_units

            M(closed + jj, closed + jj + 1) = k_on
            M(open + jj, open + jj + 1) = k_on
           #M(inactivated + i, inactivated + i + 1) = k_on

        # off rates:
            M(closed + jj + 1, closed + jj) = jj*K_C
            M(open + jj + 1, open + jj) = jj*K_O
        #  M(inactivated + i + 1, inactivated + i) = i*K_I
               

    # Transpose since we have used above: M(from, to) = k_{from, to}
    # TODO: need to check if the following is correct 
        M = np.transpose(M)
        for kk in range(n):
            M(kk,kk) = -np.sum(M, axis=0)
        
        
    #  Solve the system for Steady state at Ca_i_ss and V     
    # TODO: need to verify the following implementation   
        eq = null_space(M) #find nullspace for BK (equilibrium)
        eq = eq/np.sum(eq); #Find unique equilibrium by scaling probabilities
        #ode = @(t,y) ode_system(t,y,V)
        #[T,S] = ode15s(ode,tspan,eq)
        open = np.sum(eq[6:10]) #Calculate total steady-state open probability
        I_BK = gBK*open*(V-E_K) #Calculate steady-state current in pA
        I_K_Ca_act = I_BK
        
    else:
        I_K_Ca_act = 0.0
    
    return I_K_Ca_act
    
def potassiumPump(V, K_i, K_o):
    enable_I_K_ATP = params.enable_I_K_ATP
    if (enable_I_K_ATP == True):
        sigma   = 0.6
        g_0 = 0.05; #Testing, 3/16/16
        #g_0     = 30.95/40; % FIXME: Somewhat arbitrary. Scaled this down to match Zhou/Ferrero.
        p_0     = 0.91
        H_K_ATP = -0.001
        K_m_ATP = 0.56
        surf    = 1 # not used, what is this variable ?

        V_0 = params.V_0
        ADP_i = 10
        ATP_i = V - V_0 + ADP_i; # FIXME: arbitrary

        H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i)
        K_m = 35.8 + 17.9*ADP_i**(K_m_ATP)
        f_ATP = 1.0/(1.0 + (ATP_i/K_m)**H)

        z_K = params.z_K
        E_K = nernstPotential(z_K, K_i, K_o)
        I_K_ATP = sigma*g_0*p_0*f_ATP*(V - E_K)
    else:
        I_K_ATP = 0.0
    
    return I_K_ATP


# External stimulation
def externalStimulation(t):
    enable_I_stim = params.enable_I_stim
    if (enable_I_stim == True):
        t_cycle = params.t_cycle, t_stim = params.t_stim, I_stim_bar = params.I_stim_bar
        I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle)
    else :
        I_stim = 0.0
    
    return I_stim

# FIXME: Implement the voltage-activated hydrogen channel
def voltageActivatedHydrogen():
    enable_I_ASIC = params.enable_I_ASIC
    if (enable_I_ASIC == True):
        I_ASIC = 0.0
    else:
        I_ASIC = 0.0
    
    return I_ASIC

# Implement the TRPV4 channel
def TripCurrent(V):
    enable_I_TRPV4 = params.enable_I_TRPV4
    if (enable_I_TRPV4 == True):
        g_TRPV4 = params.g_TRPV4,  a_TRPV4 = params.a_TRPV4, b_TRPV4 = params.b_TRPV4
        if(V < 0):
            I_TRPV4 = g_TRPV4*(b_TRPV4*V + (1 - b_TRPV4)*a_TRPV4*(1 - (1 - (V/a_TRPV4))*(1 - (V/a_TRPV4))*(1-(V/a_TRPV4))))
        else:
            I_TRPV4 = 2*g_TRPV4*V**3
    else:
        I_TRPV4 = 0.0

    return I_TRPV4

# FIXME: Implement the stretch-activated TRP channel    
def stretchActivatedTrip(V):
    enable_I_TRP1 = params.enable_I_TRP1
    if (enable_I_TRP1 == True):
        g_TRP1 = params.g_TRP1, a_TRP1 = params.a_TRP1, b_TRP1 = params.b_TRP1
        if(V < 0):
            I_TRP1 = g_TRP1*(b_TRP1*V + (1 - b_TRP1)*a_TRP1*(1 - (1 - (V/a_TRP1))*(1 - (V/a_TRP1))*(1-(V/a_TRP1))))
        else:
            I_TRP1 = 2*g_TRP1*V**3
    else:
        I_TRP1 = 0.0

    return I_TRP1

# FIXME: Implement the osteo-arthritic TRP channel
def osteoArthriticTrip():
    enable_I_TRP2 = params.enable_I_TRP2
    if (enable_I_TRP2 == True):
        I_TRP2 = 0.0
    else:
        I_TRP2 = 0.0
