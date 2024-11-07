clc;

% Inputs
material = "graphite_epoxy_1";
[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", material);
schedule = [0, 0, 20, -20, 0, 90];
schedule_full = [schedule, flip(schedule)];
symmetrical = true;

if symmetrical
    n_layers = length(schedule) * 2;
else
    n_layers = length(schedule);
end

z_c = 0.005; % m

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

% Print Material and Strength Parameters (removed for brevity)

fprintf("\n================= GEOMETRY PARAMETERS =================\n\n");

fprintf('Layer Number    Type         Thickness (%s)   Orientation (degrees)\n', "mm");
fprintf('-----------------------------------------------------------------\n');

% Initialize layer number
layer_num = 1;

% Iterate through layup schedule (forward)
for i = 1:length(schedule)
    fprintf('%-13d   %-10s   %-15.3f   %20d\n', layer_num, 'ply', h_o.value * 1e3, schedule(i));
    layer_num = layer_num + 1;
end

% If symmetrical, print the core thickness
if symmetrical
    fprintf('%-13s   %-10s   %-15.3f   %20s\n', '-', 'core', 2 * z_c * 1e3, 'N/A');
    
    % Continue with the next ply layers, reverse the loop for symmetrical layup
    for i = length(schedule):-1:1
        fprintf('%-13d   %-10s   %-15.3f   %20d\n', layer_num, 'ply', h_o.value * 1e3, schedule(i));
        layer_num = layer_num + 1;
    end
end

fprintf("\n================= MATRICES =================\n");

A_matrix = calculateAMatrix(schedule, h_o.value, Q);
a_matrix = inv(A_matrix);

fprintf('\n[A] Matrix (in %s):\n', "N/m");
fprintf('%15.3e %15.3e %15.3e\n', A_matrix(1, 1), A_matrix(1, 2), A_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', A_matrix(2, 1), A_matrix(2, 2), A_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', A_matrix(3, 1), A_matrix(3, 2), A_matrix(3, 3));

fprintf('\n[a] Matrix (in %s):\n', "m/N");
fprintf('%15.3e %15.3e %15.3e\n', a_matrix(1, 1), a_matrix(1, 2), a_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', a_matrix(2, 1), a_matrix(2, 2), a_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', a_matrix(3, 1), a_matrix(3, 2), a_matrix(3, 3));

% Skateboard three-point bend inputs
P = -200 * 9.8; % [N]
L = 0.5; % [m]
b = 0.1; % [m]

M_1 = (P * L) / (4 * b);

D_matrix = calculateDMatrix(schedule, h_o.value, Q, z_c);
d_matrix = inv(D_matrix);

fprintf('\n[D] Matrix (in %s):\n', "Nm");
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(1, 1), D_matrix(1, 2), D_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(2, 1), D_matrix(2, 2), D_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(3, 1), D_matrix(3, 2), D_matrix(3, 3));

fprintf('\n[d] Matrix (in %s):\n', "(Nm)^-1");
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(1, 1), d_matrix(1, 2), d_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(2, 1), d_matrix(2, 2), d_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(3, 1), d_matrix(3, 2), d_matrix(3, 3));

fprintf("\n================= CURVATURES AND OFF-AXIS STRAIN =================\n\n");

% Skateboard three-point bend inputs
P = -200 * 9.8; % [N]
L = 0.5; % [m]
b = 0.1; % [m]
M_1 = (P * L) / (4 * b);

% D Matrix and d_matrix calculation (removed for brevity)

N_vector = [0; 0; 0];
M_vector = [M_1; 0; 0];

% Calculate curvatures
k_1 = d_matrix(1, 1) * M_vector(1);
k_2 = d_matrix(2, 1) * M_vector(1);
k_6 = d_matrix(3, 1) * M_vector(1);

k_vector = [k_1; k_2; k_6];

% Compute off-axis strain
epsilon_o_vector = a_matrix * N_vector;

% Present the variables neatly
fprintf('N_vector (N) = [%0.3f; %0.3f; %0.3f]\n\n', N_vector);
fprintf('M_vector (Nm) = [%0.3f; %0.3f; %0.3f]\n\n', M_vector);
fprintf('epsilon_o_vector = [\n  %0.3e\n  %0.3e\n  %0.3e\n]\n', epsilon_o_vector);
fprintf('\nk_vector (m^-1) = [\n  %0.3e\n  %0.3e\n  %0.3e\n]\n', k_vector);

fprintf("\n================= PER-LAYER STRESSES AND STRAINS =================\n\n");

% Function to compute on-axis and off-axis strains and stresses and return them for printing
function [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, schedule_angle, Q)
    % Calculate on-axis stress and strain vectors for bottom and top
    on_axis_strain_bottom = transformStrain(epsilon_vector_bottom, schedule_angle);
    on_axis_stress_bottom = Q * on_axis_strain_bottom;

    on_axis_strain_top = transformStrain(epsilon_vector_top, schedule_angle);
    on_axis_stress_top = Q * on_axis_strain_top;
end

% Initialize h_i and total plies
h_i = h_o.value; % Thickness of each ply
num_plies = length(schedule); % Number of plies in one side (before symmetry)

% Start at the top ply, which is located at z_c + total height of the schedule
z_i = z_c + num_plies * h_i; % Top of first ply

% Store the epsilon x values for max comparison
on_axis_epsilon_x_values = []; % Stores epsilon x values
ply_info = {}; % Stores which ply and whether it's the top or bottom

% Preallocate table for neat output
output_table = [];
stress_values = []; % To store stress values for plotting
z_heights = [];     % To store z-height values for plotting

% Loop through the first half of the layup (top side of symmetry)
for i = 1:num_plies
    angle = schedule(i); % Ply angle
    
    % z_i is the top of the ply, z_i1 is the bottom of the ply
    z_i1 = z_i - h_i; % Bottom of ply

    % Calculate strain vectors for the top and bottom of the ply
    epsilon_vector_top = epsilon_o_vector + z_i * k_vector;
    epsilon_vector_bottom = epsilon_o_vector + z_i1 * k_vector;

    % Compute strain and stress for top and bottom of ply
    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    % Append results for both Top and Bottom to the output table with a single z_height
    output_table = [output_table; {i, angle, z_i, 'Top', on_axis_strain_top(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_top)}];
    output_table = [output_table; {i, angle, z_i1, 'Bottom', on_axis_strain_bottom(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_bottom)}];

    % Store values for plotting
    stress_values = [stress_values, on_axis_stress_top(1), on_axis_stress_bottom(1)];
    z_heights = [z_heights, z_i, z_i1];

    % Store epsilon x (first element of on-axis strain) for bottom and top of ply
    on_axis_epsilon_x_values = [on_axis_epsilon_x_values, on_axis_strain_bottom(1)];
    ply_info{end+1} = sprintf('Ply %d Bottom', i);

    on_axis_epsilon_x_values = [on_axis_epsilon_x_values, on_axis_strain_top(1)];
    ply_info{end+1} = sprintf('Ply %d Top', i);

    % Update z_i for the next iteration (moving down)
    z_i = z_i1;
end

% Jump by 2 * z_c after reaching the core
z_i = z_i - 2 * z_c; % Jump by the core thickness

% Loop through the second half of the layup (bottom side of symmetry), reverse ply order
for i = num_plies:-1:1
    angle = schedule(i); % Ply angle (schedule is symmetric)

    % z_i1 is the bottom of the ply, z_i is the top
    z_i1 = z_i - h_i; % Bottom of ply

    % Calculate strain vectors for the top and bottom of the ply
    epsilon_vector_top = epsilon_o_vector + z_i * k_vector;
    epsilon_vector_bottom = epsilon_o_vector + z_i1 * k_vector;

    % Compute strain and stress for top and bottom of ply
    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    % Append results for both Top and Bottom to the output table with a single z_height
    output_table = [output_table; {num_plies + (num_plies - i + 1), angle, z_i, 'Top', on_axis_strain_top(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_top)}];
    output_table = [output_table; {num_plies + (num_plies - i + 1), angle, z_i1, 'Bottom', on_axis_strain_bottom(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_bottom)}];

    % Store values for plotting
    stress_values = [stress_values, on_axis_stress_top(1), on_axis_stress_bottom(1)];
    z_heights = [z_heights, z_i, z_i1];

    % Store epsilon x (first element of on-axis strain) for bottom and top of ply
    on_axis_epsilon_x_values = [on_axis_epsilon_x_values, on_axis_strain_bottom(1)];
    ply_info{end+1} = sprintf('Ply %d Bottom', num_plies + (num_plies - i + 1));

    on_axis_epsilon_x_values = [on_axis_epsilon_x_values, on_axis_strain_top(1)];
    ply_info{end+1} = sprintf('Ply %d Top', num_plies + (num_plies - i + 1));

    % Update z_i for the next iteration (moving down)
    z_i = z_i1;
end

% Find the overall highest magnitude of epsilon x across all plies
[max_epsilon_x, max_index] = max(abs(on_axis_epsilon_x_values));

% Convert table results to a table with units
output_table = cell2table(output_table, 'VariableNames', {'Ply', 'Angle (deg)', 'z_height (m)', 'Surface', 'Epsilon_x', 'epsilon_vector', 'on_axis_strain', 'on_axis_stress'});

% Convert the Surface column to a categorical array to avoid quotes
output_table.Surface = categorical(output_table.Surface);

% Display the table neatly
disp(output_table);

% Plotting Stress vs. Z-Height
figure;
plot(stress_values, z_heights, 'LineWidth', 2); % Remove markers and set thicker line
xlabel('Stress (MPa)');
ylabel('Z-Height (m)');
title('Stress vs Z-Height');
grid on;

fprintf("\n================= SKATEBOARD RESULTS =================\n\n");

% Display P, L, b, and calculations
fprintf('P = %.2f N\n', P);
fprintf('L = %.2f m\n', L);
fprintf('b = %.2f m\n', b);
fprintf('M_1 = (P * L) / (4 * b) = %.2f Nm\n', M_1);
fprintf('d_11 = %.3e Nm^-1\n', d_matrix(1,1));

delta_midpoint = ((P * L^3) / (48 * b)) * d_matrix(1, 1);
fprintf('delta_midpoint = ((P * L^3) / (48 * b)) * d_11 = %.5f cm\n', delta_midpoint * 1e2);

% Display the maximum epsilon x (in magnitude) and the corresponding ply and surface
fprintf("\nThe highest magnitude of epsilon_x is %f at %s.\n\n", abs(on_axis_epsilon_x_values(max_index)), ply_info{max_index});

% Display delta_midpoint right after
fprintf('The midpoint deflection (delta_midpoint) is %.5f cm.\n', delta_midpoint * 1e2);
