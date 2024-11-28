clc;

% Define constants
L = 0.2;
W = 0.1;

% Skateboard three-point bend inputs
P_1 = 1000; % [N]
L_ = 0.5; % [m]
b = 0.1; % [m]
M_1 = -(P_1 * L_) / (4 * b);
N_1 = (0.5 * P_1) / b;

M_i = [M_1; 0; 0];
N_i = [N_1; 0; 0];

% Inputs
material_string = "graphite_epoxy_1";
schedule = [0, 0, 20, -20, 0, 90];

z_c = 0.005; % m

% Evaluate the laminate
[maxstress_Rmin, quad_Rmin, hashin_Rmin, mass] = ...
    evaluateLaminate(material_string, schedule, z_c, M_i, N_i, L, W);

% Display the results cleanly
fprintf('Results of evaluateLaminate:\n');
fprintf('----------------------------------\n');
fprintf('Maximum Stress Reserve Minimum: %.4f\n', maxstress_Rmin);
fprintf('Quadratic Reserve Minimum: %.4f\n', quad_Rmin);
fprintf('Hashin Reserve Minimum: %.4f\n', hashin_Rmin);
fprintf('Mass of Laminate: %.4f g\n', mass);