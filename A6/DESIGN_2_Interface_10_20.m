clc;
clear;

% Inputs
material = "flaxpreg";
[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", material);

% Ply angle choices
ply_angles = [0, 5, -5, 20, -20, 30, -30, 50, -50, 60, -60, 90];

% Initialize results
results = [];
batch_size = 10000; % Save progress after every 10,000 iterations
batch_results = [];

% Sweep through 10 to 20 plies (symmetric, even numbers only)
for num_plies = 10:2:100
    half_plies = num_plies / 2; % Half plies for symmetric laminate
    
    % Loop over all combinations of 5 angles out of the 12 choices
    angle_combinations = nchoosek(ply_angles, 5);
    total_iterations = size(angle_combinations, 1) * factorial(5);
    current_iteration = 0;
    
    for i = 1:size(angle_combinations, 1)
        chosen_angles = angle_combinations(i, :);
        
        % Generate all permutations of the 5 chosen angles for the current number of plies
        permuted_angles = perms(chosen_angles);
        for perm_idx = 1:size(permuted_angles, 1)
            chosen_schedule = permuted_angles(perm_idx, :);
            
            % Construct symmetric schedule
            if length(chosen_schedule) <= half_plies
                % Repeat chosen schedule and truncate if needed to fit half_plies
                repeated_schedule = repmat(chosen_schedule, 1, ceil(half_plies / length(chosen_schedule)));
                half_schedule = repeated_schedule(1:half_plies); % Trim to exact half-plies
                schedule = [half_schedule, fliplr(half_schedule)]; % Symmetric laminate
            else
                continue; % Skip invalid schedules
            end
            
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
            
            % Store intermediate results with precision
            batch_results = [batch_results; {schedule, num_plies, D11, D22, round(K, 7), abs(K - 4)}];
            
            % Update progress
            current_iteration = current_iteration + 1;
            
            % Save results after every batch_size iterations
            if mod(current_iteration, batch_size) == 0
                results = [results; batch_results];
                save('intermediate_results.mat', 'results');
                batch_results = []; % Clear batch results
            end
        end
    end
end

% Append any remaining results in the final batch
if ~isempty(batch_results)
    results = [results; batch_results];
    save('intermediate_results.mat', 'results');
end

% Sort results by closeness to K = 4
results = sortrows(results, 6); % Sort by absolute difference from 4

% Display top 10 results
fprintf("\nTop 10 Results Closest to K = 4:\n");
fprintf("Schedule\t\tPlies\tD11\t\t\tD22\t\t\tK (D11/D22)\n");
for i = 1:min(10, size(results, 1))
    schedule = results{i, 1};
    num_plies = results{i, 2};
    D11 = results{i, 3};
    D22 = results{i, 4};
    K = results{i, 5};
    fprintf("%s\t%d\t%.3e\t%.3e\t%.7f\n", mat2str(schedule), num_plies, D11, D22, K);
end
