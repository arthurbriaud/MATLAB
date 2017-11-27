clear;
close all;
clc;

%------------------------------------------------------------------------------------------%
% code used for Rayleigh, mineral changes, etc... computations.
%------------------------------------------------------------------------------------------%

%% undimensionalized parameters form the simulation

mu_ref = 1e20;  % ReferenceV
k = 1e-6;       % thermdiff
DT = 1300 ;     % ReferenceT
H = 3000e3;     % layerd
T0 = 0;         % Tsurf
g = 9.8;        % gravacc
ro0= 3300;      % density
row = 1000; 
roC = 2700;
alfa = 3e-5;    % thermexp
visc_UM = 1e19; % Viscosity wanted at the surface 

%% Rayleigh numbers
disp('-------------------');
disp(' Rayleigh numbers');
disp('-------------------');
Ra = (ro0 * alfa * g * DT * (H^3))/(k * mu_ref);
                                                    
Rb1 = ((250) * g * (H^3))/(k * mu_ref);
Rb2 = ((350) * g * (H^3))/(k * mu_ref);

disp('Ra');
disp(Ra)
disp('Rb At 410 km depth');
disp(Rb1)
disp('Rb At 670 km depth');
disp(Rb2)
% For the continental lithosphere
Rb3 = ((roC -ro0) * g * (H^3))/(k * mu_ref);
disp('Ra continental litoshphere');
disp(Rb3);
%% Clayperon slopes
disp('-------------------');
disp(' Clapeyron Phases');
disp('-------------------');
gamma_1 = 3.5; % MPa/K
gamma_2 = -3;% MPa/K
h       = 3000;
Clap1 = (gamma_1 /(3.3 * g * h)) * DT;
Clap2 = (gamma_2 /(3.3* g * h)) * DT;
disp('At 410 km depth');
disp(Clap1)
disp('At 670 km depth');
disp(Clap2)

%% Size transition zone 
disp('-------------------');
disp(' size of the TZs');
disp('-------------------');
Size = 20/h;
disp (Size)

%% Check upper mantle viscosity 
Visc_UM = (visc_UM/mu_ref)*mu_ref;
disp('------------------------');
disp('viscosity at the surface');
disp('------------------------');
disp(Visc_UM)
