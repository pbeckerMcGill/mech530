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
P_1 = 1000; % [N]
L = 0.5; % [m]
b = 0.1; % [m]
M_1 = -(P_1 * L) / (4 * b);
N_1 = (0.5 * P_1) / b;

M_i = [M_1; 0; 0];
N_i = [N_1; 0; 0];

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


% Calculate curvatures
k_vector = d_matrix * M_i;
epsilon_o_vector = a_matrix * N_i;

% Display vectors
fprintf('N_vector (N) = [%0.3f; %0.3f; %0.3f]\n\n', N_i);
fprintf('M_vector (Nm) = [%0.3f; %0.3f; %0.3f]\n\n', M_i);
fprintf('epsilon_o_vector = [\n  %0.3e\n  %0.3e\n  %0.3e\n]\n', epsilon_o_vector);
fprintf('\nk_vector (m^-1) = [\n  %0.3e\n  %0.3e\n  %0.3e\n]\n', k_vector);

fprintf("\n================= PER-LAYER STRESSES AND STRAINS =================\n\n");

% Function to compute on-axis and off-axis strains and stresses
function [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, schedule_angle, Q)
    on_axis_strain_bottom = transformStrain(epsilon_vector_bottom, schedule_angle);
    on_axis_stress_bottom = Q * on_axis_strain_bottom;

    on_axis_strain_top = transformStrain(epsilon_vector_top, schedule_angle);
    on_axis_stress_top = Q * on_axis_strain_top;
end

% Initialize h_i and total plies
h_i = h_o.value;
num_plies = length(schedule);
z_i = z_c + num_plies * h_i;

output_table = {}; % Output table for layer stresses and strains

% Loop for the layup layers
for i = 1:num_plies
    angle = schedule(i);
    z_i1 = z_i - h_i;

    epsilon_vector_top = epsilon_o_vector + z_i * k_vector;
    epsilon_vector_bottom = epsilon_o_vector + z_i1 * k_vector;

    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    % Add data to output table for top and bottom of ply
    output_table = [output_table; {i, angle, z_i, 'Top', on_axis_strain_top(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_top)}];
    output_table = [output_table; {i, angle, z_i1, 'Bottom', on_axis_strain_bottom(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_bottom)}];

    % Move to the next ply
    z_i = z_i1;
end

% Add bottom side of symmetry (reversed schedule)
z_i = z_i - 2 * z_c;

for i = num_plies:-1:1
    angle = schedule(i);
    z_i1 = z_i - h_i;

    epsilon_vector_top = epsilon_o_vector + z_i * k_vector;
    epsilon_vector_bottom = epsilon_o_vector + z_i1 * k_vector;

    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    output_table = [output_table; {num_plies + (num_plies - i + 1), angle, z_i, 'Top', on_axis_strain_top(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_top), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_top)}];
    output_table = [output_table; {num_plies + (num_plies - i + 1), angle, z_i1, 'Bottom', on_axis_strain_bottom(1), sprintf('[%.3e, %.3e, %.3e]', epsilon_vector_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_strain_bottom), sprintf('[%.3e, %.3e, %.3e]', on_axis_stress_bottom)}];

    % Move to the next ply
    z_i = z_i1;
end

% Convert cell array to table for display
output_table = cell2table(output_table, 'VariableNames', {'Ply', 'Angle (deg)', 'z_height (m)', 'Surface', 'Epsilon_x', 'epsilon_vector', 'on_axis_strain', 'on_axis_stress'});

% Display the table
disp(output_table);

fprintf("\n================= QUADRATIC FAILURE CRITERIA COEFFICIENTS =================\n");

% Calculate failure criteria coefficients
F_xx = 1 / (X_t.value * X_c.value);
F_x = (1 / X_t.value) - (1 / X_c.value);
F_yy = 1 / (Y_t.value * Y_c.value);
F_y = (1 / Y_t.value) - (1 / Y_c.value);
F_s = 1 / S_c.value^2;
F_xy = sqrt(F_xx * F_yy) * (-1/2);

% Display failure criteria coefficients
fprintf('F_xx: %.3e\n', F_xx);
fprintf('F_x: %.3e\n', F_x);
fprintf('F_yy: %.3e\n', F_yy);
fprintf('F_y: %.3e\n', F_y);
fprintf('F_s: %.3e\n', F_s);
fprintf('F_xy: %.3e\n', F_xy);

fprintf("\n================= PER-LAYER STRESSES AND STRAINS =================\n\n");

% Initialize h_i and total plies
h_i = h_o.value;
num_plies = length(schedule);
z_i = z_c + num_plies * h_i;

output_table = [];
R_values_table_max_stress = [];
R_values_table_quadratic = [];

% Initialize tables as empty cell arrays with the correct number of columns
output_table = {}; % Empty cell array for output_table
R_values_table_max_stress = {}; % Empty cell array for max stress
R_values_table_quadratic = {}; % Empty cell array for quadratic criterion
R_values_table_hashin = {}; % Empty cell array for Hashin criterion

% Loop for top plies (first half)
for i = 1:num_plies
    angle = schedule(i);
    z_i1 = z_i - h_i;

    epsilon_vector_top = a_matrix * N_i + z_i * d_matrix * M_i;
    epsilon_vector_bottom = a_matrix * N_i + z_i1 * d_matrix * M_i;

    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    % Maximum stress criterion for the top surface
    R_1_top = X_t.value / on_axis_stress_top(1);
    R_2_top = abs(X_c.value / on_axis_stress_top(1));
    R_3_top = Y_t.value / on_axis_stress_top(2);
    R_4_top = abs(Y_c.value / on_axis_stress_top(2));
    R_5_top = abs(S_c.value / on_axis_stress_top(3));

    if R_1_top < 1, R_1_top = "N/A"; end
    if R_2_top < 1, R_2_top = "N/A"; end
    if R_3_top < 1, R_3_top = "N/A"; end
    if R_4_top < 1, R_4_top = "N/A"; end
    if R_5_top < 1, R_5_top = "N/A"; end

    R_values_table_max_stress = [R_values_table_max_stress; {i, angle, 'Top', R_1_top, R_2_top, R_3_top, R_4_top, R_5_top}];
    
    % Quadratic failure criterion for the top surface
    a = F_xx * on_axis_stress_top(1)^2 + 2 * F_xy * on_axis_stress_top(1) * on_axis_stress_top(2) + F_yy * on_axis_stress_top(2)^2 + F_s * on_axis_stress_top(3)^2;
    b = F_x * on_axis_stress_top(1) + F_y * on_axis_stress_top(2);
    c = -1;
    discriminant = b^2 - 4 * a * c;

    if discriminant >= 0
        R_top_1 = (-b + sqrt(discriminant)) / (2 * a);
        R_top_2 = (-b - sqrt(discriminant)) / (2 * a);
    else
        R_top_1 = "N/A";
        R_top_2 = "N/A";
    end

    R_values_table_quadratic = [R_values_table_quadratic; {i, angle, 'Top', R_top_1, R_top_2}];

    % Hashin failure criteria for the top surface
    % Tensile Fiber Mode
    if on_axis_stress_top(1) > 0
        R_hashin_1 = 1 / sqrt((on_axis_stress_top(1) / X_t.value)^2 + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_1 = "N/A";
    end
    
    % Fiber Compressive Mode
    if on_axis_stress_top(1) < 0
        R_hashin_2 = X_c.value / abs(on_axis_stress_top(1));
    else
        R_hashin_2 = "N/A";
    end
    
    % Tensile Matrix Mode
    if on_axis_stress_top(2) > 0
        R_hashin_3 = 1 / sqrt((on_axis_stress_top(2) / Y_t.value)^2 + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_3 = "N/A";
    end
    
    % Compressive Matrix Mode
    if on_axis_stress_top(2) < 0
        R_hashin_4 = 1 / sqrt((on_axis_stress_top(2) / (2 * S_c.value))^2 + ((Y_c.value / (2 * S_c.value))^2 - 1) * (on_axis_stress_top(2) / Y_c.value) + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_4 = "N/A";
    end

    R_values_table_hashin = [R_values_table_hashin; {i, angle, 'Top', R_hashin_1, R_hashin_2, R_hashin_3, R_hashin_4}];

    % Append stress values to output_table for top surface
    output_table = [output_table; {i, angle, z_i, 'Top', on_axis_stress_top(1), on_axis_stress_top(2), on_axis_stress_top(3)}];
    
    % Maximum stress criterion for the bottom surface
    R_1_bottom = X_t.value / on_axis_stress_bottom(1);
    R_2_bottom = abs(X_c.value / on_axis_stress_bottom(1));
    R_3_bottom = Y_t.value / on_axis_stress_bottom(2);
    R_4_bottom = abs(Y_c.value / on_axis_stress_bottom(2));
    R_5_bottom = abs(S_c.value / on_axis_stress_bottom(3));

    if R_1_bottom < 1, R_1_bottom = "N/A"; end
    if R_2_bottom < 1, R_2_bottom = "N/A"; end
    if R_3_bottom < 1, R_3_bottom = "N/A"; end
    if R_4_bottom < 1, R_4_bottom = "N/A"; end
    if R_5_bottom < 1, R_5_bottom = "N/A"; end

    R_values_table_max_stress = [R_values_table_max_stress; {i, angle, 'Bottom', R_1_bottom, R_2_bottom, R_3_bottom, R_4_bottom, R_5_bottom}];
    
    % Quadratic failure criterion for the bottom surface
    a = F_xx * on_axis_stress_bottom(1)^2 + 2 * F_xy * on_axis_stress_bottom(1) * on_axis_stress_bottom(2) + F_yy * on_axis_stress_bottom(2)^2 + F_s * on_axis_stress_bottom(3)^2;
    b = F_x * on_axis_stress_bottom(1) + F_y * on_axis_stress_bottom(2);
    c = -1;
    discriminant = b^2 - 4 * a * c;

    if discriminant >= 0
        R_bottom_1 = (-b + sqrt(discriminant)) / (2 * a);
        R_bottom_2 = (-b - sqrt(discriminant)) / (2 * a);
    else
        R_bottom_1 = "N/A";
        R_bottom_2 = "N/A";
    end

    R_values_table_quadratic = [R_values_table_quadratic; {i, angle, 'Bottom', R_bottom_1, R_bottom_2}];

    % Hashin failure criteria for the bottom surface
    % Tensile Fiber Mode
    if on_axis_stress_bottom(1) > 0
        R_hashin_1_bottom = 1 / sqrt((on_axis_stress_bottom(1) / X_t.value)^2 + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_1_bottom = "N/A";
    end

    % Fiber Compressive Mode
    if on_axis_stress_bottom(1) < 0
        R_hashin_2_bottom = X_c.value / abs(on_axis_stress_bottom(1));
    else
        R_hashin_2_bottom = "N/A";
    end

    % Tensile Matrix Mode
    if on_axis_stress_bottom(2) > 0
        R_hashin_3_bottom = 1 / sqrt((on_axis_stress_bottom(2) / Y_t.value)^2 + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_3_bottom = "N/A";
    end

    % Compressive Matrix Mode
    if on_axis_stress_bottom(2) < 0
        R_hashin_4_bottom = 1 / sqrt((on_axis_stress_bottom(2) / (2 * S_c.value))^2 + ((Y_c.value / (2 * S_c.value))^2 - 1) * (on_axis_stress_bottom(2) / Y_c.value) + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_4_bottom = "N/A";
    end

    R_values_table_hashin = [R_values_table_hashin; {i, angle, 'Bottom', R_hashin_1_bottom, R_hashin_2_bottom, R_hashin_3_bottom, R_hashin_4_bottom}];

    % Append stress values to output_table for bottom surface
    output_table = [output_table; {i, angle, z_i1, 'Bottom', on_axis_stress_bottom(1), on_axis_stress_bottom(2), on_axis_stress_bottom(3)}];
    
    z_i = z_i1;
end

% Adjust z_i for bottom half
z_i = z_i - 2 * z_c;

% Initialize ply number to continue from the last ply in the top half
ply_number = num_plies + 1;

% Loop for bottom plies (second half)
for i = num_plies:-1:1
    angle = schedule(i);
    z_i1 = z_i - h_i;

    epsilon_vector_top = a_matrix * N_i + z_i * d_matrix * M_i;
    epsilon_vector_bottom = a_matrix * N_i + z_i1 * d_matrix * M_i;

    [on_axis_strain_bottom, on_axis_stress_bottom, on_axis_strain_top, on_axis_stress_top] = computeStrainStress(epsilon_vector_bottom, epsilon_vector_top, angle, Q);

    % Maximum stress criterion for the top surface
    R_1_top = X_t.value / on_axis_stress_top(1);
    R_2_top = abs(X_c.value / on_axis_stress_top(1));
    R_3_top = Y_t.value / on_axis_stress_top(2);
    R_4_top = abs(Y_c.value / on_axis_stress_top(2));
    R_5_top = abs(S_c.value / on_axis_stress_top(3));

    if R_1_top < 1, R_1_top = "N/A"; end
    if R_2_top < 1, R_2_top = "N/A"; end
    if R_3_top < 1, R_3_top = "N/A"; end
    if R_4_top < 1, R_4_top = "N/A"; end
    if R_5_top < 1, R_5_top = "N/A"; end

    R_values_table_max_stress = [R_values_table_max_stress; {ply_number, angle, 'Top', R_1_top, R_2_top, R_3_top, R_4_top, R_5_top}];
    
    % Quadratic failure criterion for the top surface
    a = F_xx * on_axis_stress_top(1)^2 + 2 * F_xy * on_axis_stress_top(1) * on_axis_stress_top(2) + F_yy * on_axis_stress_top(2)^2 + F_s * on_axis_stress_top(3)^2;
    b = F_x * on_axis_stress_top(1) + F_y * on_axis_stress_top(2);
    c = -1;
    discriminant = b^2 - 4 * a * c;

    if discriminant >= 0
        R_top_1 = (-b + sqrt(discriminant)) / (2 * a);
        R_top_2 = (-b - sqrt(discriminant)) / (2 * a);
    else
        R_top_1 = "N/A";
        R_top_2 = "N/A";
    end

    R_values_table_quadratic = [R_values_table_quadratic; {ply_number, angle, 'Top', R_top_1, R_top_2}];

    % Hashin failure criteria for the top surface
    % Tensile Fiber Mode
    if on_axis_stress_top(1) > 0
        R_hashin_1 = 1 / sqrt((on_axis_stress_top(1) / X_t.value)^2 + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_1 = "N/A";
    end

    % Fiber Compressive Mode
    if on_axis_stress_top(1) < 0
        R_hashin_2 = X_c.value / abs(on_axis_stress_top(1));
    else
        R_hashin_2 = "N/A";
    end

    % Tensile Matrix Mode
    if on_axis_stress_top(2) > 0
        R_hashin_3 = 1 / sqrt((on_axis_stress_top(2) / Y_t.value)^2 + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_3 = "N/A";
    end

    % Compressive Matrix Mode
    if on_axis_stress_top(2) < 0
        R_hashin_4 = 1 / sqrt((on_axis_stress_top(2) / (2 * S_c.value))^2 + ((Y_c.value / (2 * S_c.value))^2 - 1) * (on_axis_stress_top(2) / Y_c.value) + (on_axis_stress_top(3) / S_c.value)^2);
    else
        R_hashin_4 = "N/A";
    end

    R_values_table_hashin = [R_values_table_hashin; {ply_number, angle, 'Top', R_hashin_1, R_hashin_2, R_hashin_3, R_hashin_4}];

    % Append stress values to output_table for top surface
    output_table = [output_table; {ply_number, angle, z_i, 'Top', on_axis_stress_top(1), on_axis_stress_top(2), on_axis_stress_top(3)}];
    
    % Maximum stress criterion for the bottom surface
    R_1_bottom = X_t.value / on_axis_stress_bottom(1);
    R_2_bottom = abs(X_c.value / on_axis_stress_bottom(1));
    R_3_bottom = Y_t.value / on_axis_stress_bottom(2);
    R_4_bottom = abs(Y_c.value / on_axis_stress_bottom(2));
    R_5_bottom = abs(S_c.value / on_axis_stress_bottom(3));

    if R_1_bottom < 1, R_1_bottom = "N/A"; end
    if R_2_bottom < 1, R_2_bottom = "N/A"; end
    if R_3_bottom < 1, R_3_bottom = "N/A"; end
    if R_4_bottom < 1, R_4_bottom = "N/A"; end
    if R_5_bottom < 1, R_5_bottom = "N/A"; end

    R_values_table_max_stress = [R_values_table_max_stress; {ply_number, angle, 'Bottom', R_1_bottom, R_2_bottom, R_3_bottom, R_4_bottom, R_5_bottom}];
    
    % Quadratic failure criterion for the bottom surface
    a = F_xx * on_axis_stress_bottom(1)^2 + 2 * F_xy * on_axis_stress_bottom(1) * on_axis_stress_bottom(2) + F_yy * on_axis_stress_bottom(2)^2 + F_s * on_axis_stress_bottom(3)^2;
    b = F_x * on_axis_stress_bottom(1) + F_y * on_axis_stress_bottom(2);
    c = -1;
    discriminant = b^2 - 4 * a * c;

    if discriminant >= 0
        R_bottom_1 = (-b + sqrt(discriminant)) / (2 * a);
        R_bottom_2 = (-b - sqrt(discriminant)) / (2 * a);
    else
        R_bottom_1 = "N/A";
        R_bottom_2 = "N/A";
    end

    R_values_table_quadratic = [R_values_table_quadratic; {ply_number, angle, 'Bottom', R_bottom_1, R_bottom_2}];

    % Hashin failure criteria for the bottom surface
    % Tensile Fiber Mode
    if on_axis_stress_bottom(1) > 0
        R_hashin_1_bottom = 1 / sqrt((on_axis_stress_bottom(1) / X_t.value)^2 + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_1_bottom = "N/A";
    end

    % Fiber Compressive Mode
    if on_axis_stress_bottom(1) < 0
        R_hashin_2_bottom = X_c.value / abs(on_axis_stress_bottom(1));
    else
        R_hashin_2_bottom = "N/A";
    end

    % Tensile Matrix Mode
    if on_axis_stress_bottom(2) > 0
        R_hashin_3_bottom = 1 / sqrt((on_axis_stress_bottom(2) / Y_t.value)^2 + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_3_bottom = "N/A";
    end

    % Compressive Matrix Mode
    if on_axis_stress_bottom(2) < 0
        R_hashin_4_bottom = 1 / sqrt((on_axis_stress_bottom(2) / (2 * S_c.value))^2 + ((Y_c.value / (2 * S_c.value))^2 - 1) * (on_axis_stress_bottom(2) / Y_c.value) + (on_axis_stress_bottom(3) / S_c.value)^2);
    else
        R_hashin_4_bottom = "N/A";
    end

    R_values_table_hashin = [R_values_table_hashin; {ply_number, angle, 'Bottom', R_hashin_1_bottom, R_hashin_2_bottom, R_hashin_3_bottom, R_hashin_4_bottom}];

    % Append stress values to output_table for bottom surface
    output_table = [output_table; {ply_number, angle, z_i1, 'Bottom', on_axis_stress_bottom(1), on_axis_stress_bottom(2), on_axis_stress_bottom(3)}];
    
    % Increment the ply number for continuous numbering
    ply_number = ply_number + 1;
    z_i = z_i1;
end


% Convert cell arrays to tables after loop completion for display
output_table = cell2table(output_table, 'VariableNames', {'Ply', 'Angle (deg)', 'z_height (m)', 'Surface', 'Sigma_x', 'Sigma_y', 'Sigma_s'});
R_values_table_max_stress = cell2table(R_values_table_max_stress, 'VariableNames', {'Ply', 'Angle (deg)', 'Surface', 'R_1', 'R_2', 'R_3', 'R_4', 'R_5'});
R_values_table_quadratic = cell2table(R_values_table_quadratic, 'VariableNames', {'Ply', 'Angle (deg)', 'Surface', 'R_quadratic_1', 'R_quadratic_2'});
R_values_table_hashin = cell2table(R_values_table_hashin, 'VariableNames', {'Ply', 'Angle (deg)', 'Surface', 'R_hashin_1', 'R_hashin_2', 'R_hashin_3', 'R_hashin_4'});

% Display tables
disp(output_table);
disp(R_values_table_max_stress);
disp(R_values_table_quadratic);
disp(R_values_table_hashin);

% ================= CONCLUSIONS =================

fprintf("\n================= MAX STRESS CRITERION CONCLUSIONS =================\n");

% Find the minimum R value and its layer for Max Stress criterion
min_R = inf;
min_row = [];
min_col = 0;
failure_modes = {'Fiber Tension', 'Fiber Compression', 'Matrix Tension', 'Matrix Compression', 'Shear'};

for col = 4:8
    if col > width(R_values_table_max_stress), break; end % Check column existence
    R_values_column = R_values_table_max_stress{:, col};
    % Ensure we treat the column as a cell array
    if ~iscell(R_values_column)
        R_values_column = num2cell(R_values_column);
    end
    numeric_R_values = cellfun(@(x) isnumeric(x) && ~isnan(x), R_values_column);
    if any(numeric_R_values)
        [current_min_R, idx] = min(cell2mat(R_values_column(numeric_R_values)));
        if current_min_R < min_R
            min_R = current_min_R;
            min_row = R_values_table_max_stress(numeric_R_values, :);
            min_row = min_row(idx, :); % Keep as table row
            min_col = col - 3;
        end
    end
end

failure_mode = failure_modes{min_col};
fprintf('First failure occurs in Ply %d at %.3f degrees on the %s surface.\n', min_row.Ply, min_row.("Angle (deg)"), char(min_row.Surface));
fprintf('Smallest R value (failure): R = %.3f\n', min_R);
fprintf('Failure Mode: %s\n', failure_mode);

% Calculate M_i * R and N_i * R
M_i_R = M_i * min_R;
N_i_R = N_i * min_R;
fprintf('M_i * R = [%0.3f; %0.3f; %0.3f] Nm\n', M_i_R);
fprintf('N_i * R = [%0.3f; %0.3f; %0.3f] N\n', N_i_R);

fprintf("\n================= QUADRATIC CRITERION CONCLUSIONS =================\n");

% Find the minimum positive R value and its layer for Quadratic criterion
min_R_quad = inf;
min_row_quad = [];
for col = 4:5
    if col > width(R_values_table_quadratic), break; end % Check column existence
    R_values_column = R_values_table_quadratic{:, col};
    if ~iscell(R_values_column)
        R_values_column = num2cell(R_values_column);
    end
    numeric_R_values = cellfun(@(x) isnumeric(x) && x > 0, R_values_column); % Only positive values
    if any(numeric_R_values)
        [current_min_R, idx] = min(cell2mat(R_values_column(numeric_R_values)));
        if current_min_R < min_R_quad
            min_R_quad = current_min_R;
            min_row_quad = R_values_table_quadratic(numeric_R_values, :);
            min_row_quad = min_row_quad(idx, :); % Keep as table row
        end
    end
end

fprintf('First failure occurs in Ply %d at %.3f degrees on the %s surface.\n', min_row_quad.Ply, min_row_quad.("Angle (deg)"), char(min_row_quad.Surface));
fprintf('Smallest R value (failure): R = %.3f\n', min_R_quad);

% Calculate M_i * R and N_i * R for Quadratic criterion
M_i_R_quad = M_i * min_R_quad;
N_i_R_quad = N_i * min_R_quad;
fprintf('M_i * R = [%0.3f; %0.3f; %0.3f] Nm\n', M_i_R_quad);
fprintf('N_i * R = [%0.3f; %0.3f; %0.3f] N\n', N_i_R_quad);

fprintf("\n================= HASHIN CRITERION CONCLUSIONS =================\n");

% Initialize variables to track the minimum real positive R value
min_R_hashin = inf;
min_row_hashin = [];
min_col_hashin = 0;
failure_modes_hashin = {'Fiber Tension', 'Fiber Compression', 'Matrix Tension', 'Matrix Compression'};

% Loop through each cell in columns 4 to 7 of R_values_table_hashin
for row = 1:height(R_values_table_hashin)
    for col = 4:7  % Columns with R_hashin values
        R_value = double(R_values_table_hashin{row, col});
        
        % Ignore any "N/A" or complex values, only consider real, positive numbers
        if isnumeric(R_value) && isreal(R_value) && R_value > 0
            if R_value < min_R_hashin
                min_R_hashin = R_value;
                min_row_hashin = R_values_table_hashin(row, :);  % Store entire row for context
                min_col_hashin = col - 3;  % Adjust to match failure_modes_hashin
            end
        end
    end
end

% Determine the failure mode based on the column index
failure_mode_hashin = failure_modes_hashin{min_col_hashin};

% Display the results
fprintf('First failure occurs in Ply %d at %.3f degrees on the %s surface.\n', min_row_hashin.Ply, min_row_hashin.("Angle (deg)"), char(min_row_hashin.Surface));
fprintf('Smallest R value (failure): R = %.3f\n', min_R_hashin);
fprintf('Failure Mode: %s\n', failure_mode_hashin);

% Calculate M_i * R and N_i * R for Hashin criterion
M_i_R_hashin = M_i * min_R_hashin;
N_i_R_hashin = N_i * min_R_hashin;
fprintf('M_i * R = [%0.3f; %0.3f; %0.3f] Nm\n', M_i_R_hashin);
fprintf('N_i * R = [%0.3f; %0.3f; %0.3f] N\n', N_i_R_hashin);
