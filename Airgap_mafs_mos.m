
%Calculations for ETD44 as 1 and ETD54 as 2 
Aw1 = (32.5-15.5)*(16.1e-6);
Aw2 = (41.2-18.9)*(20.2e-6);
Ac1 = 173e-6;
Ac2 = 280e-6;
Area_prod_1 = Ac1*Aw1;
Area_prod_2 = Ac2*Aw2;
Turn_len_1 = 77e-3;
Turn_len_2 = 96e-3;
Vol_1 = 2*22.3*45*15.2e-9;
At_1 = 2*(2*22.3*45+45*15.2+2*22.3*15.2)*1e-2;
At_2 = 2*(2*54.5*27.6+27.6*18.92+2*54.5*18.9)*1e-2;
Vol_2 = 2*54.5*27.6*18.9e-9;

% Given values
fs = 80e3;
q = 58*10^6;%conductivity of cu
D = 0.7;
L = 0.7875e-3; % Inductance in Henries (0.7875 mH)
I_max = 1.5; % Maximum current in Amperes
I_ripple = 0.5;
I_rms = 1.3583;
J_rms = 5e6;
k = 0.5; %packing factor 
B_max = 0.228; % Maximum flux density in Tesla (228 mT)
B_sat = 0.380; % Saturation flux density in Tesla (380 mT)
Ac = Ac1; % Core cross-sectional area in m² (173 mm²)
Aw = Aw1; % Window area in m² (278.53 mm²)
l_mag = 175.64e-3; % Magnetic path length in meters (175.64 mm)
mu_0 = 4 * pi * 1e-7; % Permeability of free space (H/m)
mu_r = 5000; % Relative permeability of core material

% Step 1: Ac Aw as Area product
Area_prod = L*I_max*I_rms/(k*B_max*J_rms);

% Step 2: Compute required number of turns (N)
N = L*I_max/(B_max*Ac);

% Step 3: Compute Air Gap (g) if needed

g = N*N*mu_0*Ac/(2*L);

% Display results
fprintf('Number of Turns (N): %.2f\n', N);
fprintf('Required Air Gap (g): %.6f meters (%.6f mm)\n', g, g * 1e3);

%Wire Calculations 
Wire_area = I_rms/J_rms;
Wire_rad = sqrt(Wire_area/pi)

%Loss Calculations
L_eff = Turn_len_1;
R_cond=  1/(q*pi*Wire_rad^2);
%R_cond = N*L_eff*1.72e-8;
Pcond = I_rms*I_rms*R_cond

%Core losses
k_d= 0.85;
alpha = 1.55;
beta = 2.75;
del_B = L*I_ripple/(N*Ac);
% Define the integral function
integral_func = @(theta) (abs(cos(theta)).^alpha) .* (2.^(beta - alpha));
% Compute the integral from 0 to 2*pi
integral_value = integral(integral_func, 0, 2*pi);
% Compute k_i
%ki = k_d / ((2*pi)^(alpha - 1) * integral_value);
ki = k_d;
P_core = Vol_1*ki*((del_B)^(beta-alpha))*fs*((abs(del_B*fs/D)^alpha)*D/fs + (abs(del_B*fs/(1-D))^alpha)*(1-D)/fs)

function P_loss = skineffect_80kHz
% This function calculates the skin effect, AC resistance, and power loss
% for a cylindrical copper wire at 80 kHz.

%% *Input Data for a Copper Wire*
u0          = 4 * pi * 1e-7;  % Permeability constant [Vs/Am]
ur          = 1;              % Relative permeability of material (Copper)
rho         = 1.72e-8;        % Resistivity of copper [Ohm*m]
D_wire_m    = 2.94e-4*2;      % Wire diameter [m]
f           = 80e3;           % Frequency = 80 kHz [Hz]

% Additional parameters
R_m         = D_wire_m / 2;   % Wire radius [m]
area        = pi * R_m^2;     % Cross-sectional area [m²]
sigma       = 1 / rho;        % Conductivity [S/m]

R_DC = 1 / (sigma * area);    % DC resistance per unit length [Ohm/m]
omega  = 2 * pi * f;         % Angular frequency [rad/s]
delta  = 1 / sqrt(omega * sigma * u0 * ur / 2);  % Skin depth [m]

%% *Compute Impedance of the Wire*
z_math = Z_wire(omega, R_m, sigma, u0, ur);

%% *Power Loss Calculation*
length_m    = 2.31;  % Length of wire in meters
current_A   = 1.25;  % RMS current in Amps

% Compute power loss at 80 kHz
P_loss = (current_A^2) * real(z_math) * length_m;  % Power loss in Watts

%% *Display Results*
fprintf('--- Skin Effect Analysis at 80 kHz ---\n');
fprintf('Skin Depth: %.6f mm\n', delta * 1e3);
fprintf('AC Resistance: %.6f Ohm/m\n', real(z_math));
fprintf('Inductance: %.6f uH/m\n', imag(z_math) / omega * 1e6);
fprintf('Power Loss: %.3f W\n', P_loss);

end

%% *Function to Compute Impedance Using Bessel Functions*
function z_math = Z_wire(omega, R, sigma, u0, ur)
  A = sqrt(-1j * omega * sigma * u0 * ur);  % Normalized argument for Bessel function
  J0 = besselj(0, A .* R);                  % Bessel function of zeroth order
  J1 = besselj(1, A .* R);                  % Bessel function of first order
  z_math = A ./ (sigma * 2 * pi .* R) .* J0 ./ J1;  % Specific impedance
end

P_skin=skineffect_80kHz;
t = 1/(sqrt(pi*mu_0*q*fs));
Z = 2*Wire_rad/(sqrt(2)*t);
GR = compute_GR(Z, 2*Wire_rad); 
H = ((N*(I_rms/l_mag)));
P_proximity = (Turn_len_1*N)*R_cond*GR*H^2



%Temperature increase
P_tot = Pcond+P_core+P_skin+P_proximity;
At = At_1;
del_T = (P_tot/At)^0.833


function GR = compute_GR(xi, d)
    % This function computes the value of GR based on the given formula, without a negative sign.
    % Inputs:
    %   xi - the variable 'ξ' in the formula
    %   d  - the variable 'd' in the formula
    % Output:
    %   GR - the calculated value based on the formula
    
    % Constants
    sqrt_2 = sqrt(2);
    
    % Compute ber and bei for orders 0, 1, 2
    ber0 = real(besselk(0, xi));
    ber1 = real(besselk(1, xi));
    ber2 = real(besselk(2, xi));
    
    bei0 = imag(besselk(0, xi));
    bei1 = imag(besselk(1, xi));
    bei2 = imag(besselk(2, xi));
    
    % First term in the parentheses
    term1 = (ber2 * ber1 + ber2 * bei1) / (ber0^2 + bei0^2);
    
    % Second term in the parentheses
    term2 = (bei2 * bei1 - bei2 * ber1) / (ber0^2 + bei0^2);
    
    % Combine the terms
    GR = (xi * pi^2 * d^2) / (2 * sqrt_2) * (term1 + term2);
end