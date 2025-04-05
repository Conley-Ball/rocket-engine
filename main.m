%% Regeneratively Cooled Rocket Engine
% Conley Ball
% April 4th, 2025

clc; close all; clear;

%% Input

thrust      = 2000; % lbf
% mass_flow   = 5; % kg/s
P_c         = 400; % psi(g)

% Propellant
fuel            = ['ethanol','water'];
fuel_ratio      = [0.75,0.25];
oxidizer        = 'LOX';
oxidizer_ratio  = 1;
O_F             = 1.5;

% Geometry
D_c = 4.5;
% L_c = 5;
L_star = 40;
conv_angle = 45;

h_ch =      [0.002 0.001 0.001 0.001];
w_rib =     [0.001 0.002 0.001 0.001];
t_ins =     [0.001 0.001 0.002 0.001];
t_out =     [0.002 0.002 0.002 0.002];
w_ch_min =  0.001;

% Other
P_a = 14.7;
T_a = 298;
c_star_eff = 0.9;

R = 8314.4598;

% Conversions

P_a = P_a*6894.76;
P_c = P_c*6894.76;
D_c = D_c*0.0254;

%% CEA

if exist('thrust','var')
    thrust = thrust * 4.4482;

    if exist('mass_flow','var')
        error('Thrust and mass flow rate cannot both be specified')
    end

    mass_flow = 1;
    thermo=CEA('problem','rocket','equilibrium','fac','ma,kg/s',mass_flow,'o/f',O_F,'p(psi)',(P_c+P_a)/6894.76,'pi/p',(P_c+P_a)/P_a,'reactants','fuel','C2H5OH(L)','wt%',fuel_ratio(1)*100,'t(k)',T_a,'fuel','H2O(L)','wt%',fuel_ratio(2)*100,'t(k)',T_a,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');

    c_star = thermo.output.eql.cstar(1)*c_star_eff;
    c_tau = thermo.output.eql.cf(end);
    mass_flow = thrust/(c_star*c_tau);

else

    if not(exist('mass_flow','var'))
        error('Either thrust or mass flow rate must be specified')
    end

    thermo=CEA('problem','rocket','equilibrium','fac','ma,kg/s',mass_flow,'o/f',O_F,'p(psi)',(P_c+P_a)/6894.76,'pi/p',(P_c+P_a)/P_a,'reactants','fuel','C2H5OH(L)','wt%',fuel_ratio(1)*100,'t(k)',T_a,'fuel','H2O(L)','wt%',fuel_ratio(2)*100,'t(k)',T_a,'oxid','O2(L)','wt%',100,'t(k)',90.0,'output','transport','mks','end');

    c_star = thermo.output.eql.cstar(1)*c_star_eff;
    c_tau = thermo.output.eql.cf(end);
    thrust = mass_flow*c_star*c_tau;

end

A_t = thrust/(P_c*c_tau);
A_e = A_t*thermo.output.eql.aeat(end);

mass_flow_f = 1/(1+O_F) * mass_flow;
mass_flow_ox = O_F/(1+O_F) * mass_flow;

fprintf('\n=============== CEA ===============\n')
fprintf('Chamber pressure:      %.1f psi(g)\n',P_c/6894.76)
fprintf('Thrust:                  %.1f lbf\n',thrust/4.4482)
fprintf('Mass flow rate:          %.3f kg/s\n',mass_flow)
fprintf('Fuel mass flow rate:     %.3f kg/s\n',mass_flow_f)
fprintf('Oxidizer mass flow rate: %.3f kg/s\n',mass_flow_ox)
fprintf('O/F Ratio:                     %0.2f\n',O_F)
fprintf('C star:                  %.1f m/s\n',c_star)
fprintf('C tau:                        %.3f\n',c_tau)
fprintf('Throat area:           %f m^2\n',A_t)
fprintf('Exit area:             %f m^2\n',A_e)

%% Cantera

% Chamber

MW_eth = 46.069;
MW_H2O = 18.015;
MW_ox = 31.998;

Y_ox = O_F/(1+O_F);
Y_eth = fuel_ratio(1)/(1+O_F);
Y_H2O = fuel_ratio(2)/(1+O_F);

X_ox = Y_ox/MW_ox;
X_H2O = Y_H2O/MW_H2O;
X_eth = Y_eth/MW_eth;

X_tot = X_ox + X_H2O + X_eth;

X_ox = X_ox/X_tot;
X_H2O = X_H2O/X_tot;
X_eth = X_eth/X_tot;

water = Solution('reactants.yaml','liquid_water');
set(water,'T',T_a,'P',(P_c+P_a));

ethanol = Solution('reactants.yaml','liquid_ethanol');
set(ethanol,'T',T_a,'P',(P_c+P_a));

ox = Solution('reactants.yaml','liquid_oxygen');
set(ox,'T',90,'P',(P_c+P_a));

gas = Solution('sandiego20161214.yaml');
mix = Mixture({water, X_H2O; ethanol, X_eth; ox, X_ox; gas, 0});
equilibrate(mix,'HP');

[gamma_ct, cp_ct,dlvpt_ct,dlvtp_ct] = get_thermo(gas);
gamma_CEA = thermo.output.eql.gamma(1);
cp_CEA = thermo.output.eql.cp(1)*1000;

% Throat

P_t = (P_c+P_a)/((gamma_ct+1)/2)^(gamma_ct/(gamma_ct-1));
M = 1;

gamma_t = gamma_ct;
error = 1;
while abs(error) > 1e-6

P_t = P_t*(1+gamma_t*M^2)/(1+gamma_t);

tgas = Solution('sandiego20161214.yaml');
set(tgas,'S',entropy_mass(gas),'P',P_t,'X',moleFractions(gas));
equilibrate(tgas,'SP');
u_t = sqrt(2*(enthalpy_mass(gas)-enthalpy_mass(tgas)));
a_t = sqrt(gamma_t*R/meanMolecularWeight(tgas)*temperature(tgas));
M = u_t/a_t;
[gamma_t,cp_t] = get_thermo(tgas);
error = 1-1/M^2;

end

% Exit

P_ratio = (P_c+P_a)/P_a;
P_e = P_a;

A_mdot = R*temperature(tgas)/(pressure(tgas)*u_t*meanMolecularWeight(tgas));

egas = Solution('sandiego20161214.yaml');

set(egas,'S',entropy_mass(tgas),'P',P_e,'X',moleFractions(tgas));
equilibrate(egas,'SP');

[gamma_e,cp_e] = get_thermo(egas);

u_e = sqrt(2*(enthalpy_mass(gas)-enthalpy_mass(egas)));
a_e = sqrt(gamma_e*R/meanMolecularWeight(egas)*temperature(egas));

A_ratio = R*temperature(egas)/(pressure(egas)*u_e*meanMolecularWeight(egas)) / A_mdot;

c_star_ct = (1/gamma_ct*((gamma_ct+1)/2)^((gamma_ct+1)/(gamma_ct-1))*1000*8.3145/meanMolecularWeight(gas)*temperature(gas))^0.5;

c_tau_ct = u_e/c_star_ct;

c_star_ct = c_star_ct*c_star_eff;

mass_flow_ct = thrust/(c_star_ct*c_tau_ct);

A_t_ct = A_mdot*mass_flow_ct;
A_t_ct = thrust/(P_c*c_tau_ct);

A_e_ct = A_t_ct*A_ratio;

fprintf('\n============= Cantera ==============\n')
fprintf('Chamber pressure:      %.1f psi(g)\n',P_c/6894.76)
fprintf('Thrust:                  %.1f lbf\n',thrust/4.4482)
fprintf('Mass flow rate:          %.3f kg/s\n',mass_flow_ct)
fprintf('Fuel mass flow rate:     %.3f kg/s\n',mass_flow_ct/(1+O_F))
fprintf('Oxidizer mass flow rate: %.3f kg/s\n',mass_flow_ct*O_F/(1+O_F))
fprintf('O/F Ratio:                     %0.2f\n',O_F)
fprintf('C star:                  %.1f m/s\n',c_star_ct)
fprintf('C tau:                        %.3f\n',c_tau_ct)
fprintf('Throat area:           %f m^2\n',A_t_ct)
fprintf('Exit area:             %f m^2\n',A_e_ct)

% figure(1)
% x = [1 2 3 4];
% hold on
% scatter(x,thermo.output.eql.gamma,'r','Marker','x')
% scatter(x,[gamma_ct,gamma_ct,gamma_t,gamma_e],'r','Marker','+')
% scatter(x,thermo.output.eql.cp,'b','Marker','x')
% scatter(x,[cp_ct,cp_ct,cp_t,cp_e]/1000,'b','Marker','+')
% hold off