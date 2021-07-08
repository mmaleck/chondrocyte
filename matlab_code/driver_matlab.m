% # This file is part of the chondrocyte modelling project at Simula
% # Research Laboratory, Norway. Refer to the files README and COPYING for
% # more information about the project as well as terms of distribution.
% #
% # Copyright (C) 2010--2016 M. M. Maleckar & Harish Narayanan
% # Licensed under the GNU GPL Version 3

% # Load (initial) model parameters
parameters;

% % Define initial conditions
x0 = [V_0, Na_i_0, K_i_0, Ca_i_0, H_i_0, Cl_i_0, a_ur_0, i_ur_0, vol_i_0, cal_0];
%x0 = [V_0, Na_i_0, K_i_0, Ca_i_0];


% % Define the initial values of parameters to be estimated
global theta0;
theta0 = [g_K_b_bar, P_K, Gmax];

% % Define the problem time range
t = linspace(0, t_final, t_final/dt);

%chondrocyte_matlab(0, x0);

%solve
%options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn',@tofile); 
%options = odeset('AbsTol', 1e-6, 'RelTol', 1e-3, 'OutputFcn',@odeprint); 
%[t,x]=ode15s(@chondrocyte_matlab, t, x0, options); 
[t,x]=ode15s(@chondrocyte_matlab, t, x0);

global ramp_Vm;
if (ramp_Vm == true)
     global V_0, global V_final, global t_final;
    V = V_0 + (V_final - V_0)*t/t_final;
% elseif 
%     (clamp_Vm == true)
%     global V_0, global V_final, global t_final;
%     for
else
    V = x(:, 1);
end

longth = length(V);
  
Na_i = x(:, 2);
K_i  = x(:, 3);
Ca_i = x(:, 4);
H_i  = x(:, 5);
Cl_i = x(:, 6);
a_ur = x(:, 7);
I_ur = x(:, 8);
vol_i = x(:, 9);
cal  = x(:, 10);

%test RMP plot
figure
subplot(2,2,1);
plot(t,V)
title('RMP')
xlabel('time (s)')
ylabel('Voltage (mV)')

subplot(2,2,2);
plot(t,Na_i)
title('Na_i')
xlabel('time (s)')
ylabel('Concentration (mM)')

subplot(2,2,3);
plot(t,K_i)
title('K_i')
xlabel('time (s)')
ylabel('Concentration (mM)')

subplot(2,2,4);
plot(t,Ca_i)
title('Ca_i')
xlabel('time (s)')
ylabel('Concentration (mM)')

% subplot(5,2,5);
% plot(t,H_i)
% 
% subplot(5,2,6);
% plot(t,Cl_i)
% 
% subplot(5,2,7);
% plot(t,a_ur)
% 
% subplot(5,2,8);
% plot(t,I_ur)
% 
% subplot(5,2,9);
% plot(t,vol_i)
% 
% subplot(5,2,10);
% plot(t,cal)


%get steady-state ion concentrations
Ca_i_ss = Ca_i(longth);
K_i_ss = K_i(longth);
Na_i_ss = Na_i(longth);
%Na_i_ss = 12.0;
V_RMP_ss = V(longth);

% To print currents function versus time in this hacky way

% ofile1 = ['I_K_DR.txt'];
% fid1 = fopen(ofile1,'w');
% ofile2 = ['I_NaK.txt'];
% fid2 = fopen(ofile2,'w');
% ofile3 = ['I_NaCa.txt'];
% fid3 = fopen(ofile3,'w');
% ofile4 = ['I_K_ur.txt'];
% fid4 = fopen(ofile4,'w');
% ofile5 = ['I_K_Ca_act.txt'];
% fid5 = fopen(ofile5,'w');
% ofile6 = ['I_K_ATP.txt'];
% fid6 = fopen(ofile6,'w');
% ofile7 = ['I_K_2pore.txt'];
% fid7 = fopen(ofile7,'w');

% Create outputstyle
% style = '';
% for n = 1:(0 + 0 + 1) % time, differential variables, currents and fluxes
%    style = [style ' %6.6e'];
% end
% style = [style '\n'];

% I_K_DR = zeros(longth,2); %1
% I_Na_K = zeros(longth,2); %2 
% I_NaCa = zeros(longth,2);  %3
% I_K_ur = zeros(longth,2); %4
% I_K_Ca_act = zeros(longth,2); %5
% I_K_ATP = zeros(longth,2);  %6
% I_K_2pore = zeros(longth,2); %7

% counter = 1;
% start = 1;
% for i = start:longth
%     % I_KDR vs. time
%     global enable_I_K_DR enable_I_NaK enable_I_NaCa enable_I_K_ur enable_I_K_Ca_act enable_I_K_ATP enable_I_K_2pore;
%     if (enable_I_K_DR == true)
%         global g_K_DR i_K_DR;
%         alpha_K_DR = 1.0/(1.0 + exp(-(V(i) + 26.7 + act_DR_shift)/4.1)); % Same as what I had. 3/3/16
%         %E_K        = nernstPotential(z_K, K_i, K_o);
%         E_K = -83;
%         I_K_DR = g_K_DR*alpha_K_DR*i_K_DR*(V(i) - E_K); %I added voltage dependence back in mwahahaahha
%         
%     else
%         I_K_DR = 0.0;
%     end
%     
%      % I_NaK vs. time
%     if (enable_I_NaK == true)
%         global I_NaK_bar K_NaK_K K_NaK_Na K_o;
%     I_NaK = I_NaK_bar*(K_o/(K_o + K_NaK_K)) ...
%         *power(Na_i(i),1.5)/(power(Na_i(i),1.5) + power(K_NaK_Na,1.5)) ...
%         *(V(i) + 150.0)/(V(i) + 200.0);
%     else
%         I_NaK = 0.0;
%     end
%     
%     % I_NaCa vs. time
%     if (enable_I_NaCa == true)
%         global F, global R, global T;
%         global Na_o, global Ca_o;
%         global K_NaCa, global gamma_Na, global d_NaCa;
% 
%     I_NaCa = K_NaCa*(  power(Na_i(i),3)*Ca_o*exp(gamma_Na*V(i)*F/(R*T)) ...
%                      - Na_o^3*Ca_i(i)*exp((gamma_Na - 1.0)*V(i)*F/(R*T))) ...
%              /(1.0 + d_NaCa*(Na_o^3*Ca_i(i) + power(Na_i(i),3)*Ca_o));
%     else
%         I_NaCa = 0.0;
%     end
%     
%     % I_K_ur vs. time
%     if (enable_I_K_ur == true)
%         global z_K g_K_ur K_o;
%     a_ur_inf   = 1.0/(1.0 + exp(-(V(i) + 26.7)/4.1));
%     i_ur_inf   = 1.0/(1.0 + exp((V(i) - 30.0)/10.0));
%     tau_a_ur   = 0.005/(1.0 + exp((V(i) + 5.0)/12.0));
%     tau_i_ur   = 0.59/(1.0 + exp((V(i) + 10.0)/24.0)) + 0.01;
%     E_K = (R*T)/(z_K*F)*log(K_o/K_i(i));
%     I_K_ur     = g_K_ur*a_ur(i)*(V(i) - E_K);
%     else
%         I_K_ur = 0.0;
%     end
%     
%     % I_K_Ca_act vs. time - not used; have the new Sun et al formulation in
%     % a different fileset, BK10 and solve_IBK
%      if (enable_I_K_Ca_act == true)
%         global T Zj Vhj ZL L0 KDc C D E N_channel z_K E_K_Ca_act;
%         kTe = 23.54*(T/273);
%         Lv = L0*exp((V(i)*ZL)/kTe);
%         Jv = exp(((V(i) - Vhj)*Zj)/kTe);
%         K = Ca_i/KDc;
%         P0 = (Lv*power((1+K*C+Jv*D+Jv*K*C*D*E),4))/((Lv*power(1+K*C+Jv*D+Jv*K*C*D*E,4))+(power(1+Jv+K+Jv*K*E,4)));
% %     % E_K = nernstPotential(z_K, K_i, K_o)
%         I_K_Ca_act_temp = N_channel*P0*Gmax*(V(i) - E_K_Ca_act);
%         if(I_K_Ca_act_temp < 0.0)
%          I_K_Ca_act = 0.0;
%         else
%          I_K_Ca_act = I_K_Ca_act_temp;
%         end
%      else
%          I_K_Ca_act = 0.0;
%      end
%     
%     % I_K_ATP vs. time -- need to clarify where this formulation came from
%     if (enable_I_K_ATP == true)
%         sigma   = 0.6;
%     %g_0     = 30.95/40; % FIXME: Somewhat arbitrary. Scaled this down to match Zhou/Ferrero.
%     g_0 = 0.05; %Testing, 3/16/16
%     p_0     = 0.91;
%     H_K_ATP = -0.001;
%     K_m_ATP = 0.56;
%     surf    = 1;
% 
%     global V_0;
%     ADP_i = 10;
%     ATP_i = V(i) - V_0 + ADP_i; % FIXME: Completely arbitrary
% 
%     H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i);
%     K_m = 35.8 + 17.9*power(ADP_i,K_m_ATP);
%     f_ATP = 1.0/(1.0 + power((ATP_i/K_m),H));
% 
%     global z_K;
%     global R, global T, global F;
%     E_K = -83.0;
%     %E_K = (R*T)/(z_K*F)*log(K_o/K_i(i));
%     I_K_ATP = sigma*g_0*p_0*f_ATP*(V(i) - E_K);
%     else
%         I_K_ATP = 0.0;
%     end
%      
%     % I_K_2pore vs. time
%     if (enable_I_K_2pore == true)
%        global F R T z_K I_K_2pore_0;
%     I_K_2pore = P_K*z_K^2*V(i)*F^2/(R*T)*(K_i(i) - K_o*exp(-z_K*V(i)*F/(R*T)))/(1 - exp(-z_K*V(i)*F/(R*T))) + I_K_2pore_0;
%     else
%         I_K_2pore = 0.0;
%     end
%     
%     counter = counter + 1
%     tpt = t(i);
%     %fprintf(fid1, '%4.2f %4.6f\n', tpt,I_K_DR);
% %     fprintf(fid2, '%4.2f %4.6f\n', tpt, I_NaK);
% %     fprintf(fid3, '%4.2f %4.6f\n', tpt, I_NaCa);
%      %fprintf(fid4, '%4.2f %4.6f\n', tpt, I_K_ur);
% %     fprintf(fid5, '%4.2f %4.6f\n', tpt, I_K_Ca_act);
% %     fprintf(fid6, '%4.2f %4.6f\n', tpt, I_K_ATP);
% %     fprintf(fid7, '%4.2f %4.6f\n', tpt, I_K_2pore);
% end


%%%% VOLTAGE CLAMP IV PRINT %%%
% Run a voltage clamp for currents ; SS concentrations are
% used and not evolved to avoid interpolation, for the time being

ofile8 = ['IV_K2pore.txt'];
fid8 = fopen(ofile8,'w');
ofile9 = ['voltage.txt'];
fid9 = fopen(ofile9,'w');
ofile10 = ['IV_KATP.txt'];
fid10 = fopen(ofile10,'w');
ofile11 = ['IV_INab.txt'];
fid11 = fopen(ofile11,'w');
ofile12 = ['IV_IKb.txt'];
fid12 = fopen(ofile12,'w');
ofile13 = ['IV_IClb.txt'];
fid13 = fopen(ofile13,'w');
ofile14 = ['IV_INaK.txt'];
fid14 = fopen(ofile14,'w');
ofile15 = ['IV_INaCa.txt'];
fid15 = fopen(ofile15,'w');
ofile16 = ['IV_ICaP.txt'];
fid16 = fopen(ofile16,'w');
ofile17 = ['IV_I_TRPV4.txt'];
fid17 = fopen(ofile17,'w');
ofile18 = ['IV_I_bg.txt'];
fid18 = fopen(ofile18,'w');
ofile19 = ['IV_I_total.txt'];
fid19 = fopen(ofile19,'w');
ofile20 = ['IV_I_K_Ca.txt'];
fid20 = fopen(ofile20,'w');
ofile21 = ['IV_I_K_DR.txt'];
fid21 = fopen(ofile21,'w');
ofile22 = ['IV_act_DR.txt'];
fid22 = fopen(ofile22,'w');
ofile23 = ['RMP_vs_KDR.txt'];
fid23 = fopen(ofile23,'w');
ofile24 = ['IV_TRPV4.txt'];
fid24 = fopen(ofile24,'w');
ofile25 = ['IV_RMP.txt'];
fid25 = fopen(ofile25,'w');

VV = -150:0.1:100; %voltage clamp in mV
    for i = 1:length(VV)
        global VV
        
        
        % I_K_DR (printed in pA/pF)
        global g_K_DR i_K_DR act_DR_shift;
        alpha_K_DR(i) = 1.0/(1.0 + exp(-(VV(i) + 26.7 + act_DR_shift)/4.1)); % Same as what I had. 3/3/16
        %E_K        = nernstPotential(z_K, K_i, K_o);
        E_K = -83;
        I_K_DR(i) = g_K_DR*alpha_K_DR(i)*i_K_DR*(VV(i) - E_K); %pA
        fprintf(fid21, '%4.2f %4.6f\n', VV(i), I_K_DR(i)/C_m); %pA/pF
        fprintf(fid22, '%4.2f %4.6f\n', VV(i), alpha_K_DR(i));
        
        % I_Na_K (in pA; printed IV pA/pF)
        global I_NaK_bar K_NaK_K K_NaK_Na K_o Na_i_clamp;
        %K_o = 10;
        I_NaK(i) = I_NaK_bar*(K_o/(K_o + K_NaK_K)) ...
        *power(Na_i_clamp,1.5)/(power(Na_i_clamp,1.5) + power(K_NaK_Na,1.5)) ...
        *(VV(i) + 150.0)/(VV(i) + 200.0);
        %fprintf(fid14, '%4.2f %4.6f\n', VV(i), I_NaK(i));
        fprintf(fid14, '%4.2f %4.6f\n', VV(i), I_NaK(i)/C_m);
        
        % I_NaCa (in pA; printed IV pA/pF)
        global F, global R, global T;
        global Na_o, global Ca_o;
        global K_NaCa, global gamma_Na, global d_NaCa;
        global Na_i_0; 
        global Ca_i_0;
        global Na_i_clamp;
        global NCX_scale;
        
        I_NaCa(i) = NCX_scale*K_NaCa*(Na_i_clamp^3*Ca_o*exp(gamma_Na*VV(i)*F/(R*T)) ...
                    - Na_o^3*Ca_i_0*exp((gamma_Na - 1.0)*VV(i)*F/(R*T))) ...
            /(1.0 + d_NaCa*(Na_o^3*Ca_i_0 +Na_i_clamp^3*Ca_o));
        
        % Below uses scalars for IC variables to check and print I_NCX
        % function
        %I_NaCa(i) = NCX_scale*K_NaCa*(  12.0^3*Ca_o*exp(gamma_Na*VV(i)*F/(R*T)) ...
         %            - Na_o^3*1.8*exp((gamma_Na - 1.0)*VV(i)*F/(R*T))) ...
         %    /(1.0 + d_NaCa*(Na_o^3*1.8 +12.0^3*Ca_o));
        
%         %I_NaCa(i) = K_NaCa*((power(Na_i_0,3)*Ca_o*exp(gamma_Na*VV(i)*F/(R*T)) ...
%                      - power(Na_o,3)*Ca_i_0*exp((gamma_Na - 1.0)*VV(i)*F/(R*T))) ...
%              /(1.0 + d_NaCa*(power(Na_o,3)*Ca_i_0 + power(Na_i_0,3)*Ca_o)));
         %fprintf(fid15, '%4.2f %4.6f\n', VV(i), I_NaCa(i));
         fprintf(fid15, '%4.2f %4.6f\n', VV(i), I_NaCa(i)/C_m);
         
         % I_Ca_ATP (pA)
         global I_Ca_ATP_bar, global k_Ca_ATP;
         I_Ca_ATP(i) = I_Ca_ATP_bar*(Ca_i_ss/(Ca_i_ss + k_Ca_ATP)); %using SS value of Ca_i from ODE solve above
         fprintf(fid16, '%4.2f %4.6f\n', VV(i), I_Ca_ATP(i));
         
         % I_K_ATP (pA?) Zhou/Ferrero, Biophys J, 2009
        sigma   = 0.6;
        g_0     = 30.95/40; % FIXME: Somewhat arbitrary. Scaled this down to match Zhou/Ferrero.
        %g_0 = 30.95; %Testing, 3/16/16
        p_0     = 0.91;
        
        H_K_ATP = -0.001;
        K_m_ATP = 0.56;
        %surf    = 10^4; Zhou, for cardiomyocyte
        surf = 1;
        
        global V_0;
        ADP_i = 10;
        ATP_i = VV(i) - V_0 + ADP_i; % FIXME: Completely arbitrary
        %ATP_i = 1000;
            
        H = 1.3 + 0.74*exp(-H_K_ATP*ADP_i);
        K_m = 35.8 + 17.9*power(ADP_i,K_m_ATP);
        f_ATP = 1.0/(1.0 + power((ATP_i/K_m),H));

        global z_K;
        global R, global T, global F;
        %E_K = (R*T)/(z_K*F)*log(K_o/K_i_0)
        E_K = -94.02; % From Zhou, et al
        I_K_ATP(i) = sigma*g_0*p_0*f_ATP*(VV(i) - E_K);
        fprintf(fid10,'%4.2f %4.6f\n', VV(i), I_K_ATP(i));
        
        % I_K_2pore modeled as a simple Boltzmann
        % relationship via GHK, scaled to match isotonic K+ data from Bob Clark (pA; pA/pF in print)
        %K_i_2pore = K_i_ss;
        %K_i_2pore = K_i_0;
        global F R T z_K I_K_2pore_0 Q K_i_0;
        %I_K_2pore(i) = 1.0*(P_K)*z_K^2*VV(i)*F^2/(R*T)*(K_i_2pore - K_o*exp(-z_K*VV(i)*F/(R*T)))/(1 - exp(-z_K*VV(i)*F/(R*T))) + I_K_2pore_0;
        I_K_2pore(i) = 5.0*Q*sqrt(K_o/K_i_0)*VV(i)*(1 - (K_o/K_i_0)*exp(-z_K*VV(i)*F/(R*T)))/(1- exp(-z_K*VV(i)*F/(R*T))) + I_K_2pore_0;

        fprintf(fid8, '%4.2f %4.6f\n', VV(i), I_K_2pore(i)/C_m);
        fprintf(fid9, '%4.2f\n', VV(i));
        
        % I_Na_b (pA; pA/pF in print)
        global z_Na, global g_Na_b_bar, global Na_o_0;
        %E_Na = nernstPotential(z_Na, Na_i_0, Na_o_0);
        E_Na = 55.0;
        I_Na_b(i) = g_Na_b_bar*(VV(i) - E_Na);
        fprintf(fid11, '%4.2f %4.6f\n', VV(i), I_Na_b(i)/C_m);
        
        % I_K_b (pA; pA/pF in print)
        global z_K;
        %E_K = nernstPotential(z_K, K_i, K_o);
        E_K = -83;
        I_K_b(i) = g_K_b_bar*(VV(i) - E_K);
        fprintf(fid12, '%4.2f %4.6f\n', VV(i), I_K_b(i)/C_m);
        
        % I_Cl_b (pA; pA/pF in print)
        global z_Cl, global g_Cl_b_bar, global Cl_o, global C_m;
        %E_Cl = nernstPotential(z_Cl, Cl_o, Cl_i);
        E_Cl = -65.0; % -40.0 was previously used?
        I_Cl_b(i) = g_Cl_b_bar*(VV(i) - E_Cl);
        %I_Cl_b(i) = 0.0;
        fprintf(fid13, '%4.2f %4.6f\n', VV(i), I_Cl_b(i)/C_m);
        
        % I_leak (pA); not printed, added to I_bg
        global g_leak;
        I_leak(i) = g_leak*VV(i);
    
        
        % I_TRPV4 (pA)
        global g_TRPV4, global a_TRPV4, global b_TRPV4;
        
        I_TRPV4(i) = 0.0;
        %I_TRPV4(i) = 1/50*g_TRPV4*VV(i)^3;
%         I_TRPV4(i) = g_TRPV4*(VV(i));
%         if(VV(i) < 0)
%             I_TRPV4(i) = g_TRPV4*(b_TRPV4*VV(i) + (1 - b_TRPV4)*a_TRPV4*(1 - (1 - (VV(i)/a_TRPV4))*(1 - (VV(i)/a_TRPV4))*(1-(VV(i)/a_TRPV4))));
%         else
%             I_TRPV4(i) = 2*g_TRPV4*VV(i)^3;
%         end
        fprintf(fid17, '%4.2f %4.6f\n', VV(i), I_TRPV4(i));
        
        % I_bg (pA; pA/pF in print)
%         I_bg(i)= I_Na_b(i) + I_K_b(i) + I_Cl_b(i) + I_K_DR(i);
        I_bg(i)= (I_Na_b(i) + I_K_b(i) + I_Cl_b(i) + I_leak(i));
        fprintf(fid18, '%4.2f %4.6f\n', VV(i), I_bg(i)/C_m);
        
        
  %I_K_Ca_act (new version) (pA), with converted Ca_i units for model,
  %print as pA/pF
        % Set constants
        convert_units = 1e6; % Convert from nM (e-9) to mM (e-3)
        %Ca_i_BK = 3e-6; % Assume units of mM
        %gBK = 100.0;
        %gBK = 2.50;%Gmax, nS
        E_K = -83; %Sun, et al
        %tspan = 1e-3*[0:1:300]; %time in seconds
        %C_m = 7; %pF

        K_C = 17;
        K_O = 0.5;

        % 1/25 factor for A, nothing for B was default
        A = 1/65*[0.659,3.955,25.05,129.2,261.1];
        B = 4.5*[2651.7, 1767.8, 1244.0, 713.0, 160.0];

        z_CO = 0.718;  
        z_OC = 0.646;

        FVonRT = F*VV(i)/(R*T);

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
            k_on = (5-jj)*Ca_i_ss*convert_units;

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
        open(i) = sum(eq(6:10)); %Calculate total steady-state open probability
        I_BK(i) = gBK*open(i)*(VV(i)-E_K); %Calculate steady-state current in pA
        
        fprintf(fid20, '%4.2f %4.6f\n', VV(i), I_BK(i)/C_m);
      
        
        % I TRPV4 (pA; pA/pF in print)
        global g_TRPV4, global a_TRPV4, global b_TRPV4;
        if(VV(i) < 0)
        I_TRPV4(i) = g_TRPV4*(b_TRPV4*VV(i) + (1 - b_TRPV4)*a_TRPV4*(1 - (1 - (VV(i)/a_TRPV4))*(1 - (VV(i)/a_TRPV4))*(1-(VV(i)/a_TRPV4))));
        else
        I_TRPV4(i) = 2*g_TRPV4*VV(i)^3;
        end
        fprintf(fid24, '%4.2f %4.6f\n', VV(i), I_TRPV4(i)/C_m);
  
        
        % I_RMP (pA; pA/pF in print)
        I_RMP(i)= (I_bg(i) + I_BK(i) + I_K_DR(i) + I_NaCa(i) + I_NaK(i) + I_K_2pore(i));
        fprintf(fid25, '%4.2f %4.6f\n', VV(i), I_RMP(i)/C_m);
        
        % I_total (print in pA/pF - check all the currents summed are in pA)
        I_total(i) = (I_NaK(i) + I_NaCa(i) + I_Ca_ATP(i) ...
      + I_K_DR(i) + I_K_2pore(i) + I_K_ATP(i) + I_BK(i) ...
      + I_Na_b(i) + I_K_b(i) + I_Cl_b(i) + I_leak(i) + I_TRPV4(i))/C_m;
     
        %+ I_total
      fprintf(fid19, '%4.2f %4.6f\n', VV(i), I_total(i));
      
      
      
        
    end
 
fclose('all');
slope_G = (I_bg(length(VV))-I_bg(1))/(VV(length(VV))-VV(1)) % pA/mV = nS
R = 1/slope_G % = GOhms
Ca_i_ss;
%
