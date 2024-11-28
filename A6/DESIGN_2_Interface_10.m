clc;

% Inputs
material = "flaxpreg";
[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", material);

% Ply angle choices
ply_angles = [0, 5, -5, 20, -20, 30, -30, 50, -50, 60, -60, 90];

% Initialize results
results = [];

% Loop over all combinations of 5 angles out of the 12 choices
angle_combinations = nchoosek(ply_angles, 5);

total_iterations = size(angle_combinations, 1) * factorial(5);
current_iteration = 0;

for i = 1:size(angle_combinations, 1)
    chosen_angles = angle_combinations(i, :);
    
    % Generate all permutations of the 5 chosen angles for a 10-layer symmetric laminate
    permuted_angles = perms(chosen_angles);
    for perm_idx = 1:size(permuted_angles, 1)
        chosen_schedule = permuted_angles(perm_idx, :);
        schedule = [chosen_schedule, fliplr(chosen_schedule)];
        
        % Build full schedule
        schedule_full = [schedule, flip(schedule)];
        symmetrical = true;
        
        % Compute z_c (mid-plane shift)
        z_c = 0; % Assuming mid-plane is zero
        
        % Compute Q's
        Q_xx = m.value * E_x.value;
        Q_yy = m.value * E_y.value;
        Q_yx = m.value * nu_x.value * E_y.value;
        Q_xy = m.value * nu_y.value * E_x.value;
        Q_ss = E_s.value;
        
        Q = [Q_xx Q_xy 0; Q_yx Q_yy 0; 0 0 Q_ss];
        
        % Calculate D matrix
        D_matrix = calculateDMatrix(schedule, h_o.value, Q, z_c);
        
        % Extract D11 and D22
        D11 = D_matrix(1, 1);
        D22 = D_matrix(2, 2);
        
        % Calculate K = D11 / D22
        K = D11 / D22;
        
        % Store result
        results = [results; {schedule, D11, D22, K, abs(K - 4)}];
        
        % Update progress
        current_iteration = current_iteration + 1;
        progress_fraction = current_iteration / total_iterations;
        
        % Print progress
        fprintf("Progress: %.2f%% (%d/%d iterations)\n", progress_fraction * 100, current_iteration, total_iterations);
    end
end

% Sort results by closeness to K = 4
results = sortrows(results, 5); % Sort by absolute difference from 4

% Display top 10 results
fprintf("\nTop 10 Results Closest to K = 4:\n");
fprintf("Schedule\t\tD11\t\t\tD22\t\t\tK (D11/D22)\n");
for i = 1:min(10, size(results, 1))
    schedule = results{i, 1};
    D11 = results{i, 2};
    D22 = results{i, 3};
    K = results{i, 4};
    fprintf("%s\t%.3e\t%.3e\t%.7f\n", mat2str(schedule), D11, D22, K);
end
