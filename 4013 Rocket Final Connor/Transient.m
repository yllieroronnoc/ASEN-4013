close all
clear all
clc

t(1) = 0; % [s] initial time
rb(1) = 0; % [m] initial burn grain displacement

%% Digitized Data
T_ven = csvread('Thrust.csv');
T_ven(:,2) = T_ven(:,2);
%Plot Digitized Vendor Data
figure('Name','Digitized Data');
plot(T_ven(:,1),T_ven(:,2),'LineWidth',2);
hold on 
xlabel('Time in Seconds');
ylabel('Thrust in Pounds');
title('Vendor Provided Thrust Curve for F26T with Blue Thunder Propellant');
grid on ;
xlim([0 1.2])
ylim([0 16])


%% INPUTS
cstar_eff = .75; % [-], cstar efficiency
t_step = .03; % [s] time step
P_atm = 101325; % [Pa] ambient pressure
a = 0.000069947; % [-] burn rate coefficient
n = 0.321; % [-] burn rate exponent
cstar = 1500; % [m/s] characteristic velocity
h_grain = 1.505; % [in] motor grain height
r_grain_inner = 0.177/2; % [in] motor grain radius
r_grain_outer = 0.908/2;
r_throat = 0.123/2; % [in] throat radius
r_exit = 0.231/2; % [in] exit radius
Mass = 0.025 ; % [kg] Propellant mass

%% CONVERSIONS
h_grain = h_grain*0.0254; % [m] motor grain height
r_grain_inner = r_grain_inner*0.0254; % [m] motor grain inner radius
r_grain_outer = r_grain_outer*0.0254; % [m] motor grain outer radius
r_throat = r_throat*0.0254; % [m] throat radius
r_exit = r_exit*0.0254; % [m] exit radius

%% QUANTITY CALCULATIONS
Vol = h_grain * ( (r_grain_outer)^2 - (r_grain_inner)^2 ) * pi(); % [m^3]
rho_p = Mass/Vol; % [kg/m^3]
A_throat = pi()*(r_throat)^2; % [m^2]
A_exit = pi()*(r_exit)^2; % [m^2]
AR_sup = A_exit/A_throat; % supersonic area ratio

V_burn = 0; % [m^3]
j = 1;
while rb < r_grain_outer && rb < h_grain % while there is unburned grain remaining
    [A_burn(j)] = burn_geometry(r_grain_inner,r_grain_outer, h_grain, rb); % [m] burn area, burn cavity volume
    Pc(j) = ((a * rho_p * A_burn(j) * cstar) / (A_throat)).^((1)/(1-n))/1e6; % [MPa] chamber pressure
    burn_rate(j) = a*(Pc(j)*10^6)^n; % [m/s] burn rate
    rb = rb + burn_rate(j) * t_step; % [m] updates burn displacement
    burn_rate(j)
    % delta_Vol = ; % [m^3/s] rate of change in burn cavity volume 
    [T_predicted(j),cstar] = thrust_calc(P_atm, Pc(j), A_exit, rho_p, burn_rate(j), A_burn(j), AR_sup);
    cstar = cstar*cstar_eff; % [m/s]
    %action time 
    if j == 1
        t(j) = t_step;
    else
        t(j) = t(j-1) + t_step;
    end
    j = j+1;
end
%Plot Digitized Vendor Data
figure('Name','Digitized Data');
plot(T_ven(:,1),T_ven(:,2),'LineWidth',2);
hold on 
plot(t,T_predicted,'LineWidth',2)
xlabel('Time in Seconds');
ylabel('Thrust in Pounds');
xlim([0 1.4])
title('Adjusted Model vs Vendor Data');
grid on ;


%% Performance computations

%maximum thrust calculation
max_t_ven = max(T_ven(:,2));
max_t_mod = max(T_predicted);
%total impulse from digitized data
I_ven_tot = trapz(T_ven(:,1),T_ven(:,2));
%Convert to SI units (lb-s to N-s)
I_ven_tot = I_ven_tot * 4.44822162;
%total impulse from model 
I_mod_tot = trapz(t,T_predicted);
%Convert to SI units (lb-s to N-s)
I_mod_tot = I_mod_tot * 4.44822162;

%avg spec imp
isp_av_ven = I_ven_tot/(Mass * 9.81);
isp_av_mod = I_mod_tot/(Mass * 9.81);
%avg effective velocity
avg_cstar_ven = isp_av_ven*9.81;
avg_cstar_mod = isp_av_mod*9.81;
%action time
% 10% max T
T_10_ven = max_t_ven * .1;
%ven
for i = 1:length(T_ven)
    if(T_ven(i,2) <= T_10_ven)
       time_fin = T_ven(i,1);
       break
    end
end
act_ven = time_fin;

% 10% max T
T_10_mod = max_t_mod * .1;
k = 0;
%ven
for i = 1:length(T_predicted)
    if((k == 0) && (T_predicted(i) >= T_10_mod))
        time_int = t(i);
        k = k + 1;
    end
    if((k == 1)&& (T_predicted(i) <= T_10_mod))
       time_fin = t(i);
       break
    end
end
act_mod = time_fin - time_int;



