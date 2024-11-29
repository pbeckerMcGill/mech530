clc;

% Define constants
L = 0.2;
W = 0.1;

% Skateboard three-point bend inputs

% M_i = [-991;-99;-105];
% N_i = [-22500;3100;-2000];

M_i = [-950;-95;-115];
N_i = [-20800;-2900;-2250];
% Inputs

material_string = "graphite_epoxy_1";
schedule = [0, 30, -30, 45, -45, 30, -30, 45, -45, 90, 90];

z_c = 0.00075; % m

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