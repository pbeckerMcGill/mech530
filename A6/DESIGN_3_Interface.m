% Define constants
L = 0.2;
W = 0.1;
z_c = 0.00075;

% Loading cases
loading_cases = struct('M_i', {[-991; -99; -105], [-950; -95; -115]}, ...
                       'N_i', {[-22500; 3100; -2000], [-20800; -2900; -2250]});

% Material strings
materials = {'graphite_epoxy_1', 'fiberglass', 'kevlar_epoxy', ...
             'graphite_epoxy_2', 'flaxpreg', 'graphite_thermoplastic'};

% Ply angles with 5° increments
angles = 0:5:90; % Generate angles from 0 to 90 in 5° steps
ply_angles = [-angles(2:end-1), angles(2:end-1), 0, 90]; % Include positive and negative angles except for 0 and 90

% Maximum plies per side (symmetrical laminate)
max_half_plies = 10;

% Initialize results storage
results = [];
iteration_count = 0; % Track iterations for backups

% Total iterations (approximation for progress calculation)
num_materials = length(materials);
num_loading_cases = length(loading_cases);
num_schedules = sum(arrayfun(@(x) size(nchoosek(ply_angles, x), 1), 1:max_half_plies));
total_iterations = num_materials * num_schedules * num_loading_cases;

% Generate all schedules
schedules = {};
for num_plies = 1:max_half_plies
    perms = nchoosek(ply_angles, num_plies); % Generate combinations of ply angles
    for i = 1:size(perms, 1)
        schedules{end+1} = perms(i, :); %#ok<SAGROW>
    end
end

disp("hi")

% Iterate through materials
for material_idx = 1:length(materials)
    material_string = materials{material_idx};
    
    % Iterate through schedules
    for schedule_idx = 1:length(schedules)
        first_half_schedule = schedules{schedule_idx}; % Use only the first half of the laminate
        
        % Iterate through loading cases
        for case_idx = 1:num_loading_cases
            M_i = loading_cases(case_idx).M_i;
            N_i = loading_cases(case_idx).N_i;
            
            % Evaluate the laminate
            [maxstress_Rmin, quad_Rmin, hashin_Rmin, mass] = ...
                evaluateLaminate(material_string, first_half_schedule, z_c, M_i, N_i, L, W);
            
            % Store results
            results = [results; struct('material', material_string, ...
                                       'schedule', first_half_schedule, ...
                                       'M_i', M_i, ...
                                       'N_i', N_i, ...
                                       'maxstress_Rmin', maxstress_Rmin, ...
                                       'quad_Rmin', quad_Rmin, ...
                                       'hashin_Rmin', hashin_Rmin, ...
                                       'mass', mass)]; %#ok<SAGROW>
            
            % Increment iteration count
            iteration_count = iteration_count + 1;
            
            % Display progress
            progress_percentage = (iteration_count / total_iterations) * 100;
            fprintf('Progress: %d/%d iterations (%.2f%%)\n', iteration_count, total_iterations, progress_percentage);
            
            % Save backup every 50,000 iterations
            if mod(iteration_count, 50000) == 0
                backup_filename = sprintf('laminate_results_backup_%d.mat', iteration_count);
                save(backup_filename, 'results');
                fprintf('Backup saved: %s (Iteration %d)\n', backup_filename, iteration_count);
            end
        end
    end
end

% Final save of results
save('laminate_results.mat', 'results');
disp('Sweeps complete. Results saved to laminate_results.mat.');
