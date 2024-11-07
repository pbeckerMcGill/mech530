% Inputs
%   - Modulus Parameters
%   - Strength Parameters
%   - Geometry Parameters

material = "graphite_epoxy_2";

[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", material);
schedule = [45, -45, 20, -20, 0, 90];
symmetrical = true;

if symmetrical
    n_layers = length(schedule) * 2;
else
    n_layers = length(schedule);
end

z_c = 0.075; % m

% Outputs
%   - Input Modulus Parameters
%   - Input Strength Parameters
%   - Input Geometry Parameters
%   - On-axis S & Q Matrices

% Compute S's
S_xx = 1 / E_x.value;
S_xy = -nu_y.value / E_y.value;
S_yx = -nu_x.value / E_x.value;
S_yy = 1 / E_y.value;
S_ss = 1 / E_s.value;

% Compute Q's
Q_xx = m.value * E_x.value;
Q_yy = m.value * E_y.value;
Q_yx = m.value * nu_x.value * E_y.value;
Q_xy = m.value * nu_y.value * E_x.value;
Q_ss = E_s.value;

% Build matrices
S = [S_xx S_xy 0; S_yx S_yy 0; 0 0 S_ss];
Q = [Q_xx Q_xy 0; Q_yx Q_yy 0; 0 0 Q_ss];

fprintf('CHOSEN MATERIAL: %s\n\n', material);

% Cleaned up titles for better formatting and readability
fprintf("\n================= MATERIAL PROPERTIES =================\n\n");

fprintf("Modulus Parameters\n");
fprintf('%-5s: %15.3f %s\n', 'E_x', E_x.value, E_x.unit);
fprintf('%-5s: %15.3f %s\n', 'E_y', E_y.value, E_y.unit);
fprintf('%-5s: %15.3f %s\n', 'E_s', E_s.value, E_s.unit);
fprintf('%-5s: %15.3f %s\n', 'nu_x', nu_x.value, nu_x.unit);
fprintf('%-5s: %15.3f %s\n', 'nu_y', nu_y.value, nu_y.unit);
fprintf('%-5s: %15.3f %s\n', 'm', m.value, m.unit);

fprintf("\nStrength Parameters\n");
fprintf('%-5s: %15.3f %s\n', 'X_t', X_t.value, X_t.unit);
fprintf('%-5s: %15.3f %s\n', 'X_c', X_c.value, X_c.unit);
fprintf('%-5s: %15.3f %s\n', 'Y_t', Y_t.value, Y_t.unit);
fprintf('%-5s: %15.3f %s\n', 'Y_c', Y_c.value, Y_c.unit);
fprintf('%-5s: %15.3f %s\n', 'S_c', S_c.value, S_c.unit);
fprintf('%-5s: %15.3f %s\n', 'h_o', h_o.value * 1e3, "mm");
fprintf('%-5s: %15.3f %s\n', 'rho', rho.value, rho.unit);

fprintf("\n================= GEOMETRY PARAMETERS =================\n\n");

% Print the header for the layup table
fprintf('Layer Number   Layer Type   Thickness (%s)   Orientation (degrees)\n', "mm");
fprintf('-----------------------------------------------------------------\n');

% Initialize layer number
layer_num = 1;

% Iterate through layup schedule (forward)
for i = 1:length(schedule)
    fprintf('%-13d   %-10s   %-15.3f   %20d\n', layer_num, 'ply', h_o.value * 1e3, schedule(i));
    layer_num = layer_num + 1;
end

% If symmetrical, print the core thickness in the same table
if symmetrical
    fprintf('%-13d   %-10s   %-15.3f   %20s\n', layer_num, 'core', 2 * z_c * 1e3, 'N/A');
    layer_num = layer_num + 1;

    % Iterate through layup schedule (backward)
    for i = length(schedule):-1:1
        fprintf('%-13d   %-10s   %-15.3f   %20d\n', layer_num, 'ply', h_o.value * 1e3, schedule(i));
        layer_num = layer_num + 1;
    end
end

fprintf("\n================= ON-AXIS MATRICES =================\n");

% Directly print matrix [S] without a loop
fprintf('\nMatrix [S] (in %s^-1):\n', E_x.unit);
fprintf('%15.3e %15.3e %15.3e\n', S(1, 1), S(1, 2), S(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', S(2, 1), S(2, 2), S(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', S(3, 1), S(3, 2), S(3, 3));

% Directly print matrix [Q] without a loop
fprintf('\nMatrix [Q] (in %s):\n', E_x.unit);
fprintf('%15.3f %15.3f %15.3f\n', Q(1, 1), Q(1, 2), Q(1, 3));
fprintf('%15.3f %15.3f %15.3f\n', Q(2, 1), Q(2, 2), Q(2, 3));
fprintf('%15.3f %15.3f %15.3f\n', Q(3, 1), Q(3, 2), Q(3, 3));

fprintf("\n================= Assignment 2, Question 1 =================\n");

% Add transformation and off-axis matrices for each layer
fprintf("\n================= OFF-AXIS MATRICES PER LAYER =================\n");

for i = 1:length(schedule)
    theta = schedule(i);
    
    % Transform the Q and S matrices
    Q_transf = transformQ(Q, theta);
    S_transf = transformS(S, theta);
    
    fprintf('\nLayer %d - Orientation: %d degrees\n', i, theta);
    
    % Display the transformed [S] matrix
    fprintf('\nTransformed Matrix [S] (in %s^-1):\n', E_x.unit);
    fprintf('%15.3e %15.3e %15.3e\n', S_transf(1, 1), S_transf(1, 2), S_transf(1, 3));
    fprintf('%15.3e %15.3e %15.3e\n', S_transf(2, 1), S_transf(2, 2), S_transf(2, 3));
    fprintf('%15.3e %15.3e %15.3e\n', S_transf(3, 1), S_transf(3, 2), S_transf(3, 3));

    % Display the transformed [Q] matrix
    fprintf('\nTransformed Matrix [Q] (in %s):\n', E_x.unit);
    fprintf('%15.3f %15.3f %15.3f\n', Q_transf(1, 1), Q_transf(1, 2), Q_transf(1, 3));
    fprintf('%15.3f %15.3f %15.3f\n', Q_transf(2, 1), Q_transf(2, 2), Q_transf(2, 3));
    fprintf('%15.3f %15.3f %15.3f\n', Q_transf(3, 1), Q_transf(3, 2), Q_transf(3, 3));
end

fprintf("\n================= Assignment 2, Question 2 =================\n");

% Given stress vector for Assignment 2, Question 2
ply_angle_a2 = 30;  % Example ply angle
off_axis_stress_vector = [9990; -3100; -4400];

% Print the given off-axis stress vector
fprintf("\nOff-Axis Stress Vector (Given) (MPa):\n");
fprintf('%15.3e\n', off_axis_stress_vector);

% Transform S matrix using ply_angle_a2
S_transf_a2 = transformS(S, ply_angle_a2);

% Compute off-axis strain vector
off_axis_strain_vector = S_transf_a2 * off_axis_stress_vector;

% Print the results for off-axis strain vector
fprintf("\nOff-Axis Strain Vector:\n");
fprintf('%15.3e\n', off_axis_strain_vector);

% Calculate on-axis stress and strain vectors
on_axis_stress_vector = transformStress(off_axis_stress_vector, ply_angle_a2);
on_axis_strain_vector = transformStrain(off_axis_strain_vector, ply_angle_a2);

% Print the results for on-axis stress and strain vectors
fprintf("\nOn-Axis Stress Vector (MPa):\n");
fprintf('%15.3e\n', on_axis_stress_vector);

fprintf("\nOn-Axis Strain Vector:\n");
fprintf('%15.3e\n', on_axis_strain_vector);
