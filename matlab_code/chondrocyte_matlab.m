% # This file is part of the chondrocyte modelling project at Simula
% # Research Laboratory, Norway. Refer to the files README and COPYING for
% # more information about the project as well as terms of distribution.
% #
% # Copyright (C) 2010--2016  M. M. Maleckar & Harish Narayanan
% # Licensed under the GNU GPL Version 3
% 
% % The solution fields and parameters have the following units:
% %
% % Voltage: mV
% % Current: pA
% % Time: s
% % Concentration: mM/l
% % Conductance: pS
% % Capacitance: pF



% Define the ODE model for the electrophysiology of the chondrocyte
function xdot = ode_rhs(t, x)

  global theta0;
  xdot = ode_rhs_parametrized(x, t, theta0);

%endfunction

function xdot = ode_rhs_parametrized(x, t, theta)

  % Initialize and populate vector of unknowns
  global apply_Vm;
  if (apply_Vm == true)
    V = appliedVoltage(t);
  else
  V = x(1);
  end
  
  Na_i = x(2);
  K_i  = x(3);
  Ca_i = x(4);
  H_i  = x(5);
  Cl_i = x(6);
  a_ur = x(7);
  i_ur = x(8);
  global vol_i_0;
  vol_i = vol_i_0; %x(9); %FIXME: Fixing the volume while debugging
  cal   = x(10); %calmodulin

  % Define external concentrations
  K_o = appliedPotassiumConcentration(t);

  % Extract parameters
  g_K_b_bar = theta(1);
  P_K = theta(2);
  Gmax = theta(3);

  % Calculate background currents
  I_Na_b = backgroundSodium(V, Na_i);
  I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar);
  I_Cl_b = backgroundChloride(V, Cl_i);
  I_leak = backgroundLeak(V);

  % Calculate pump and exchanger currents
  I_NaK = sodiumPotassiumPump(V, Na_i, K_i, K_o);
  I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i);
  I_NaH = sodiumHydrogenExchanger(Na_i, H_i);
  I_Ca_ATP = calciumPump(Ca_i);

  % Calculate potassium currents
  I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur);
  I_K_DR = DelayedRectifierPotassium(V);
  I_K_2pore = twoPorePotassium(V, K_i, K_o, P_K);
  I_K_Ca_act = calciumActivatedPotassium(V, Ca_i);
  I_K_ATP = potassiumPump(V, K_i, K_o);
  
  
  % Calculate other currents
  I_ASIC = voltageActivatedHydrogen();
  I_TRP1 = stretchActivatedTrip(V);
  I_TRP2 = osteoArthriticTrip();
  I_TRPV4 = TripCurrent(V);
  I_stim = externalStimulation(t);

  % Total ionic contribution (pA)
  I_i = I_Na_b + I_K_b + I_Cl_b + I_leak ...
      + I_NaK + I_NaCa + I_Ca_ATP ...
      + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP ...
      + I_ASIC + I_TRP1 + I_TRP2 + I_TRPV4;

% I_i = I_Na_b + I_K_b + I_Cl_b ...
%         + I_NaCa + I_Ca_ATP ...
%       + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP ...
%       + I_ASIC + I_TRP1 + I_TRP2;

  % Determine incremental changes in evolving quantities
  global F;
  global C_m;

  % Evolve calmodulin - Ca buffer
  cal_dot = 200000.0*Ca_i*(1.0 - cal) - 476.0*cal;
   % cal_dot = 200000.0*Ca_i*(1.0 - cal);
   
  % Evolve the concentrations
    global clamp_Na_i;
    global clamp_K_i;
    if (clamp_Na_i == true)
        Na_i_dot = 0;
    else
        Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH + I_Cl_b + 0.5*I_leak)/(vol_i*F);
       %Na_i_dot = - (I_Na_b + 3*I_NaK + 3*I_NaCa - I_NaH)/(vol_i*F);
    end
    if (clamp_K_i == true)
        K_i_dot = 0;
    else
        K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP + 0.5*I_leak)/(vol_i*F);
        %K_i_dot  = - (I_K_b  - 2*I_NaK + I_K_ur + I_K_DR + I_K_2pore + I_K_Ca_act + I_K_ATP)/(vol_i*F);
    end
  
  %Ca_i_dot =   -(I_Ca_ATP - 2*I_NaCa + I_TRPV4)/(2*vol_i*F) - 0.045*cal_dot;
  Ca_i_dot =   -(I_Ca_ATP - 2*I_NaCa + I_TRPV4)/(2*vol_i*F) - 0.045*cal_dot;
  %H_i_dot =  - (I_NaH)/(vol_i*F);
  H_i_dot = 0;
  Cl_i_dot =  (I_Cl_b)/(vol_i*F);
  %Cl_i_dot = 0;
  
  [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);

  a_ur_dot = (a_ur_inf - a_ur)/tau_a_ur;
  i_ur_dot = (i_ur_inf - i_ur)/tau_i_ur;

  global Na_o;
  global Ca_o;
  global H_o;
  global Cl_o;

  global Na_i_0;
  global K_i_0;
  global Ca_i_0;
  global H_i_0;
  global Cl_i_0;
  
%   global Na_i_clamp;

  % Think this is volume from one of the UK papers...check later 3/16/2016
  osm_i_0 = Na_i_0 + K_i_0 + Ca_i_0 + H_i_0 + Cl_i_0;
  osm_i = Na_i + K_i + Ca_i + H_i + Cl_i;
  osm_o = Na_o + K_o + Ca_o + H_o + Cl_o;
  dosm = osm_i_0 - osm_o;

  P_f = 10.2e-4;
  SA = 6.0^(2.0/3.0)*pi^(1.0/3.0)*vol_i^(2.0/3.0);
  V_W = 18.0;
  vol_i_dot = P_f*SA*V_W*(osm_i - osm_o - dosm);

  
 %Initialize ODEs for state variables
  xdot = zeros(10, 1);

  global apply_Vm;
  if (apply_Vm == true)
    xdot(1) = 0.0;
  else
    xdot(1) = 1/C_m*(-I_i + I_stim);
  end

  global clamp_conc;
  if (clamp_conc == true)
    xdot(2) = 0.0;
    xdot(3) = 0.0;
    xdot(4) = 0.0;
    xdot(5) = 0.0;
    xdot(6) = 0.0;
  else
    xdot(2) = Na_i_dot;
    xdot(3) = K_i_dot;
    xdot(4) = Ca_i_dot;
    xdot(5) = H_i_dot;
    xdot(6) = Cl_i_dot;
  end

  xdot(7) = a_ur_dot;
  xdot(8) = i_ur_dot;
  xdot(9) = vol_i_dot;
  xdot(10) = cal_dot;

%endfunction


% Potential of an ion X across the membrane (mV).
function E_X = nernstPotential(z, X_i, X_o)
  global R, global T, global F;
  E_X = (R*T)/(z*F)*log(X_o/X_i);
%endfunction

% Applied voltage (mV)
function V = appliedVoltage(t)
  global clamp_Vm, global ramp_Vm, global step_Vm;
  if (clamp_Vm == true)
    global V_0;
    V = V_0;
  elseif (ramp_Vm == true)
    global V_0, global V_final, global t_final;
    V = V_0 + (V_final - V_0)*t/t_final;
  elseif (step_Vm == true)
    global t_cycle t_stim;
    V = (ceil((t - 30)/t_cycle).*square((t - 30)*2*pi/t_cycle, t_stim/t_cycle) + ceil((t - 30)/t_cycle))/2*10 - 90;
    if (V == 0) V = 0.01; end
  end
%endfunction

function K_o = appliedPotassiumConcentration(t)
  global step_K_o K_o_0;
  if (step_K_o == false)
    K_o = K_o_0;
  else
    if (t <= 10)
      K_o = 5;
    elseif (t > 10 && t <= 20)
      K_o = 30;
    elseif (t > 20 && t <= 30)
      K_o = 75;
    elseif (t > 30 && t <= 40)
      K_o = 140;
    else
      K_o = 5;
    end
  end
%endfunction


% Background sodium current from "Ionic channels of excitable
% membranes," B. Hille. (pA)
function I_Na_b = backgroundSodium(V, Na_i)
  global enable_I_Na_b;
  if (enable_I_Na_b == true)
    global z_Na, global g_Na_b_bar, global Na_o; 
    E_Na = nernstPotential(z_Na, Na_i, Na_o);
    I_Na_b = g_Na_b_bar*(V - E_Na);
  else
    I_Na_b = 0.0;
  end
%endfunction

% Background potassium current from "Ionic channels of excitable
% membranes," B. Hille. (pA)
function I_K_b = backgroundPotassium(V, K_i, K_o, g_K_b_bar)
  global enable_I_K_b;
  if (enable_I_K_b == true)
    global z_K;
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_b = g_K_b_bar*(V - E_K);
  else
    I_K_b = 0.0;
  end
%endfunction

% Background chloride current from "Ionic channels of excitable
% membranes," B. Hille. (pA)
function I_Cl_b = backgroundChloride(V, Cl_i)
  global enable_I_Cl_b;
  if (enable_I_Cl_b == true)
    global z_Cl, global g_Cl_b_bar, global Cl_o;
    %E_Cl = nernstPotential(z_Cl, Cl_o, Cl_i);
    %E_Cl = -40.0;
    E_Cl = -65.0;
    I_Cl_b = g_Cl_b_bar*(V - E_Cl);
  else
   I_Cl_b = 0.0;
  end
%endfunction

function I_leak = backgroundLeak(V)
  global enable_I_leak;
  if (enable_I_leak == true)
      global g_leak;
    I_leak = g_leak*V;
  else
    I_leak = 0.0;
  end


% % Sodium-potassium pump from "Mathematical Model of an Adult Human
% % Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
% % Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
% % Giles. Circ. Res. 1998; 82; 63-81 (Table 12, pp. 77) (pA)

function I_NaK = sodiumPotassiumPump(V, Na_i, K_i, K_o)
  global enable_I_NaK Na_i_0;
  if (enable_I_NaK == true)
    global I_NaK_bar K_NaK_K K_NaK_Na;
    I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) ...
        *(Na_i_0^1.5/(Na_i_0^1.5 + K_NaK_Na^1.5)) ...
        *(V + 150.0)/(V + 200.0);
  else
    I_NaK = 0.0;
  end
%endfunction

% % Sodium-calcium exchanger from "Mathematical Model of an Adult Human
% % Atrial Cell: The Role of K+ Currents in Repolarization," A. Nygren, C.
% % Fiset, L. Firek, J. W. Clark, D. S. Lindblad, R. B. Clark and W. R.
% % Giles. Circ. Res. 1998; 82; 63-81 (Table 13, pp. 77) (pA)

function I_NaCa = sodiumCalciumExchanger(V, Na_i, Ca_i)
  global enable_I_NaCa;
  if (enable_I_NaCa == true)
    global F, global R, global T;
    global Na_o, global Ca_o;
    global K_NaCa, global gamma_Na, global d_NaCa, global Na_i_0;
    global NCX_scale;

    I_NaCa = NCX_scale*K_NaCa*(  Na_i_0^3*Ca_o*exp(gamma_Na*V*F/(R*T)) ...
                     - Na_o^3*Ca_i*exp((gamma_Na - 1.0)*V*F/(R*T))) ...
             /(1.0 + d_NaCa*(Na_o^3*Ca_i + Na_i_0^3*Ca_o));
  else
    I_NaCa = 0.0;
  end
%endfunction

% % Sodium-hydrogen exchanger from "A Model of Na+/H+ Exchanger and Its
% % Central Role in Regulation of pH and Na+ in Cardiac Myocytes," Chae
% % Young Cha, Chiaki Oka, Yung E. Earm, Shigeo Wakabayashi, and Akinori
% % Noma. Biophysical Journal 2009; 97; 2674-2683 (pp. 2675) (pA)

function I_NaH = sodiumHydrogenExchanger(Na_i, H_i)
  global enable_I_NaH;
  if (enable_I_NaH == true)
    global n_H, global K_H_i_mod;
    global k1_p, global k1_m, global k2_p, global k2_m;
    global Na_o, global H_o, global N_NaH_channel;
    global K_Na_o, global K_H_o, global K_Na_i, global K_H_i;

    I_NaH_mod  = 1/(1 + (K_H_i_mod^n_H/H_i^n_H));
    t1 = k1_p*Na_o/K_Na_o / (1 + Na_o/K_Na_o + H_o/K_H_o);
    t2 = k2_p*H_i/K_H_i   / (1 + Na_i/K_Na_i + H_i/K_H_i);
    t3 = k1_m*Na_i/K_Na_i / (1 + Na_i/K_Na_i + H_i/K_H_i);
    t4 = k2_m*H_o/K_H_o   / (1 + Na_o/K_Na_o + H_o/K_H_o);
    I_NaH_exch = (t1*t2 - t3*t4) / (t1 + t2 + t3 + t4);
    I_NaH = N_NaH_channel*I_NaH_mod*I_NaH_exch;
  else
    I_NaH = 0.0;
  end
%endfunction

% %  Calcium pump from Nygren et al. (pA)

function I_Ca_ATP = calciumPump(Ca_i)
  global enable_I_Ca_ATP;
  if (enable_I_Ca_ATP == true)
     global I_Ca_ATP_bar, global k_Ca_ATP;
     I_Ca_ATP = I_Ca_ATP_bar*(Ca_i/(Ca_i + k_Ca_ATP));
  else
    I_Ca_ATP = 0.0;
  end
%endfunction


% % Ultra-rapidly rectifying potassium channel from "Action potential rate
% % dependence in the human atrial myocyte," M. M. Maleckar, J. L.
% % Greenstein, W. R. Giles and N. A. Trayanova. Am. J. Physiol. Heart.
% % Circ. Physiol. 2009; 297; 1398-1410 (Appendix, pp. 1408) (pA)

function [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V)

  a_ur_inf   = 1.0/(1.0 + exp(-(V + 26.7)/4.1));
  i_ur_inf   = 1.0/(1.0 + exp((V - 30.0)/10.0));
  tau_a_ur   = 0.005/(1.0 + exp((V + 5.0)/12.0));
  tau_i_ur   = 0.59/(1.0 + exp((V + 10.0)/24.0)) + 0.01;
%endfunction

function I_K_ur = ultrarapidlyRectifyingPotassium(V, K_i, K_o, a_ur, i_ur)
  global enable_I_K_ur;
  if (enable_I_K_ur == true)
    global z_K g_K_ur;
    [a_ur_inf, i_ur_inf, tau_a_ur, tau_i_ur] = ultraRapidlyRectifyingPotassiumHelper(V);
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur     = g_K_ur*a_ur*(V - E_K);
  else
    I_K_ur = 0.0;
  end
%endfunction


% % Delayed rectifer potassium channel from Clark, et al (pA)
function I_K_DR = DelayedRectifierPotassium(V)
    global enable_I_K_DR;
    if (enable_I_K_DR == true)
        global g_K_DR i_K_DR act_DR_shift;
        alpha_K_DR = 1.0/(1.0 + exp(-(V + 26.7 + act_DR_shift)/4.1)); % Same as what I had. 3/3/16
        %E_K        = nernstPotential(z_K, K_i, K_o);
        E_K = -83;
        I_K_DR = g_K_DR*alpha_K_DR*i_K_DR*(V - E_K);
    else
        I_K_DR = 0.0;
    end
              
        
% % From Bob Clark et al., J. Physiol. 2011, Figure 4 - IKDR
function I_K_ur_ref = ultrarapidlyRectifyingPotassium_ref(V, K_i, K_o)
    G_K = 28.9;  % pS/pF
    V_h = -26.7; % mV
    S_h = 4.1;   % mV
    global C_m z_K;
    E_K        = nernstPotential(z_K, K_i, K_o);
    I_K_ur_ref = G_K*(V - E_K)/(1 + exp(-(V - V_h)/S_h)) * C_m/1000.0;
%endfunction

% % Two-pore potassium current - novel?; modeled as a simple Boltzmann
% relationship via GHK (pA)
function I_K_2pore = twoPorePotassium(V, K_i_0, K_o, Q)
  global enable_I_K_2pore;
  if (enable_I_K_2pore == true)
    global F R T z_K I_K_2pore_0;
    %OLD, via Harish: I_K_2pore = P_K*z_K^2*V*F^2/(R*T)*(K_i_0 - K_o*exp(-z_K*V*F/(R*T)))/(1- exp(-z_K*V*F/(R*T))) + I_K_2pore_0;
    I_K_2pore = 5*Q*sqrt(K_o/K_i_0)*V*(1 - (K_o/K_i_0)*exp(-z_K*V*F/(R*T)))/(1- exp(-z_K*V*F/(R*T))) + I_K_2pore_0;
  else
    I_K_2pore = 0.0;
  end
%endfunction

% % Calcium-activated potassium current - Sun, et al formulations (pA)
% % FIXME: Check the following carefully
function I_K_Ca_act = calciumActivatedPotassium(V, Ca_i)
  global enable_I_K_Ca_act;
  if (enable_I_K_Ca_act == true)
      global F R T;
      %I_K_Ca_act (new version) (pA), with converted Ca_i units for model
        % Set constants
        convert_units = 1e6; % Convert from nM (e-9) to mM (e-3)
        gBK = 2.50;%Gmax, nS
        E_K = -83; %Sun, et al
        %tspan = 1e-3*[0:1:300]; %time in seconds
        %C_m = 7; %pF

        K_C = 17;
        K_O = 0.5;

        A = 1/65*[0.659,3.955,25.05,129.2,261.1];
        B = 4.5*[2651.7, 1767.8, 1244.0, 713.0, 160.0];

        z_CO = 0.718;  
        z_OC = 0.646;

        FVonRT = F*V/(R*T);

    for m = 1:5
        alpha_BK(m) = A(m)*exp(+z_CO*FVonRT);
        beta_BK(m)  = B(m)*exp(-z_OC*FVonRT);
        %delta(i) = C(i)*exp(+z_OI*FVonRT);
        %gamma(i) = D(i)*exp(-z_IO*FVonRT);
    end

        % Build Markov Matrix
        n = 10; 
        M = zeros(n,n); 

        % numbering scheme, offsets:
        closed = 0;
        open = 5;
        %inactivated = 10;


        % Vertical transitions
        for k = 1:5

            M(closed + k , open + k) = alpha_BK(k);
            M(open + k , closed + k) = beta_BK(k);

        %  M(open + i, inactivated + i) = delta(i);
        %  M(inactivated + i, open + i) = gamma(i);

        end

        % Horizontal transitions
        for jj = 1:4

        % on rates:
            k_on = (5-jj)*Ca_i*convert_units;

            M(closed + jj, closed + jj + 1) = k_on;
            M(open + jj, open + jj + 1) = k_on;
        % M(inactivated + i, inactivated + i + 1) = k_on;

        % off rates:
            M(closed + jj + 1, closed + jj) = jj*K_C;
            M(open + jj + 1, open + jj) = jj*K_O;
        %  M(inactivated + i + 1, inactivated + i) = i*K_I;

        end

% Transpose since we have used above: M(from, to) = k_{from, to}
        M = M';
        for kk = 1:n
        M(kk,kk) = -sum(M(:,kk));
        end
        
% Solve the system for Steady state at Ca_i_ss and V        
        eq = null(M); %find nullspace for BK (equilibrium)
        eq = eq/sum(eq); %Find unique equilibrium by scaling probabilities
        %ode = @(t,y) ode_system(t,y,V);
        %[T,S] = ode15s(ode,tspan,eq);
        open = sum(eq(6:10)); %Calculate total steady-state open probability
        I_BK = gBK*open*(V-E_K); %Calculate steady-state current in pA
        I_K_Ca_act = I_BK;
        
 %%%% OLD VERSION %%%%%%     
%     global T Zj Vhj ZL L0 KDc C D E N_channel z_K E_K_Ca_act;
%     kTe = 23.54*(T/273);
%     Lv = L0*exp((V*ZL)/kTe);
%     Jv = exp(((V - Vhj)*Zj)/kTe);
%     K = Ca_i/KDc;
%     P0 = (Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)/((Lv*(1+K*C+Jv*D+Jv*K*C*D*E)^4)+((1+Jv+K+Jv*K*E)^4));
% %     % E_K = nernstPotential(z_K, K_i, K_o)
%     I_K_Ca_act_temp = N_channel*P0*Gmax*(V - E_K_Ca_act);
%     if(I_K_Ca_act_temp < 0.0)
%       I_K_Ca_act = 0.0;
%     else
%       I_K_Ca_act = I_K_Ca_act_temp;
%     end
  else
    I_K_Ca_act = 0.0;
  end
%endfunction

% % ATP-powered potassium pump from "Modeling Cardiac Action Potential
% % Shortening Driven by Oxidative Stress-Induced Mitochondrial
% % Oscillations in Guinea Pig Cardiomyocytes," L. Zhou, S. Cortassa, A-C.
% % Wei, M. A. Aon, R. L. Winslow, and B. O'Rourke. Biophys. J. 2009; 97;
% % 1843-1852 (p. 1845)
% % FIXME: Clean up the following
% % FIXME: Check the following carefully
% % (pA)???

function I_K_ATP = potassiumPump(V, K_i, K_o)
  global enable_I_K_ATP;
  if (enable_I_K_ATP == true)
    sigma   = 0.6;
    %g_0     = 30.95/40; % FIXME: Somewhat arbitrary. Scaled this down to match Zhou/Ferrero.
    g_0 = 0.05; %Testing, 3/16/16
    p_0     = 0.91;
    H_K_ATP = -0.001;
    K_m_ATP = 0.56;
    surf    = 1;

    global V_0;
    ADP_i = 10;
    ATP_i = V - V_0 + ADP_i; % FIXME: arbitrary

    H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i);
    K_m = 35.8 + 17.9*ADP_i^K_m_ATP;
    f_ATP = 1.0/(1.0 + (ATP_i/K_m)^H);

    global z_K;
    E_K = nernstPotential(z_K, K_i, K_o);
    I_K_ATP = sigma*g_0*p_0*f_ATP*(V - E_K);
  else
    I_K_ATP = 0.0;
  end
%endfunction


% # External stimulation
function I_stim = externalStimulation(t)
  global enable_I_stim;
  if (enable_I_stim == true)
    global t_cycle, global t_stim, global I_stim_bar;
    I_stim = I_stim_bar*square(t*2*pi/t_cycle, t_stim/t_cycle);
  else
    I_stim = 0.0;
  end
%endfunction

% # FIXME: Implement the voltage-activated hydrogen channel
function I_ASIC = voltageActivatedHydrogen()
  global enable_I_ASIC;
  if (enable_I_ASIC == true)
    I_ASIC = 0.0;
  else
    I_ASIC = 0.0;
  end
%endfunction

% Implement the TRPV4 channel
function I_TRPV4 = TripCurrent(V)
  global enable_I_TRPV4;
  if (enable_I_TRPV4 == true)
    global g_TRPV4, global a_TRPV4, global b_TRPV4;
    if(V < 0)
      I_TRPV4 = g_TRPV4*(b_TRPV4*V + (1 - b_TRPV4)*a_TRPV4*(1 - (1 - (V/a_TRPV4))*(1 - (V/a_TRPV4))*(1-(V/a_TRPV4))));
    else
      I_TRPV4 = 2*g_TRPV4*V^3;
    end
  else
    I_TRPV4 = 0.0;
  end
  
  
% # FIXME: Implement the stretch-activated TRP channel
function I_TRP1 = stretchActivatedTrip(V)
  global enable_I_TRP1;
  if (enable_I_TRP1 == true)
    global g_TRP1, global a_TRP1, global b_TRP1;
    if(V < 0)
      I_TRP1 = g_TRP1*(b_TRP1*V + (1 - b_TRP1)*a_TRP1*(1 - (1 - (V/a_TRP1))*(1 - (V/a_TRP1))*(1-(V/a_TRP1))));
    else
      I_TRP1 = 2*g_TRP1*V^3;
    end
  else
    I_TRP1 = 0.0;
  end
%endfunction

% # FIXME: Implement the osteo-arthritic TRP channel
function I_TRP2 = osteoArthriticTrip()
  global enable_I_TRP2;
  if (enable_I_TRP2 == true)
    I_TRP2 = 0.0;
  else
    I_TRP2 = 0.0;
  end
%endfunction
