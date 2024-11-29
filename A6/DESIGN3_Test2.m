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

M_i_1 = [-991;-99;-105];
N_i_1 = [-22500;3100;-2000];

% Inputs
material_string = "graphite_epoxy_1";
schedule = [-25];

z_c = 0.00075; % m

% Evaluate the laminate
[maxstress_Rmin, quad_Rmin, hashin_Rmin, mass] = ...
    evaluateLaminate(material_string, schedule, z_c, M_i_1, N_i_1, L, W);

% Display the results cleanly
fprintf('Results of evaluateLaminate for load 1:\n');
fprintf('----------------------------------\n');
fprintf('Maximum Stress Reserve Minimum: %.4f\n', maxstress_Rmin);
fprintf('Quadratic Reserve Minimum: %.4f\n', quad_Rmin);
fprintf('Hashin Reserve Minimum: %.4f\n', hashin_Rmin);
fprintf('Mass of Laminate: %.4f g\n', mass);

% M_i_2 = [-991;-99;-105];
% N_i_2 = [-22500;3100;-2000];
% 
% % Evaluate the laminate
% [maxstress_Rmin, quad_Rmin, hashin_Rmin, mass] = ...
%     evaluateLaminate(material_string, schedule, z_c, M_i_2, N_i_2, L, W);
% 
% % Display the results cleanly
% fprintf('Results of evaluateLaminate for load 2:\n');
% fprintf('----------------------------------\n');
% fprintf('Maximum Stress Reserve Minimum: %.4f\n', maxstress_Rmin);
% fprintf('Quadratic Reserve Minimum: %.4f\n', quad_Rmin);
% fprintf('Hashin Reserve Minimum: %.4f\n', hashin_Rmin);
% fprintf('Mass of Laminate: %.4f g\n', mass);