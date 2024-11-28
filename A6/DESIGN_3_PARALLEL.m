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

% Sweep range: 5 to 10 plies
min_half_plies = 5;
max_half_plies = 10;

% Backup settings
backup_interval = 50000; % Save backup every 50,000 iterations
backup_folder = 'backups_parallel'; % Folder to save backups
if ~exist(backup_folder, 'dir')
    mkdir(backup_folder); % Create backup folder if it doesn't exist
end

% Parallel computation
results = {}; % Initialize results
backup_counter = 0; % Track backups

% Create parallel pool
parpool; % Start a parallel pool

% Outer loop for materials
for material_idx = 1:length(materials)
    material_string = materials{material_idx};
    
    % Outer loop for number of plies in the first half (5 to 10 plies)
    for num_plies = min_half_plies:max_half_plies
        % Generate schedules for the current number of plies
        schedules = nchoosek(ply_angles, num_plies);
        
        % Preallocate local results for parallel processing
        local_results = cell(size(schedules, 1), 1);
        
        % Parallel loop for schedules
        parfor schedule_idx = 1:size(schedules, 1)
            first_half_schedule = schedules(schedule_idx, :); % Current schedule
            
            % Temporary container for valid results
            temp_results = {};
            success = true; % Track success for this schedule
            
            % Iterate through loading cases
            for case_idx = 1:length(loading_cases)
                M_i = loading_cases(case_idx).M_i;
                N_i = loading_cases(case_idx).N_i;
                
                try
                    % Call evaluateLaminate directly (no suppression)
                    [maxstress_Rmin, quad_Rmin, hashin_Rmin, mass] = evaluateLaminate(material_string, first_half_schedule, z_c, M_i, N_i, L, W);
                    
                    % Check Rmin values for validity
                    if maxstress_Rmin >= 2 && quad_Rmin >= 2 && hashin_Rmin >= 2
                        % Store result for this loading case
                        temp_results{end+1} = struct('material', material_string, ...
                                                     'schedule', first_half_schedule, ...
                                                     'M_i', M_i, ...
                                                     'N_i', N_i, ...
                                                     'maxstress_Rmin', maxstress_Rmin, ...
                                                     'quad_Rmin', quad_Rmin, ...
                                                     'hashin_Rmin', hashin_Rmin, ...
                                                     'mass', mass);
                    else
                        success = false; % Skip if any Rmin is below 2
                        break; % Exit the loop
                    end
                catch
                    success = false; % Skip on any error
                    break; % Exit the loop
                end
            end
            
            % If successful, save results
            if success && ~isempty(temp_results)
                local_results{schedule_idx} = vertcat(temp_results{:});
            else
                local_results{schedule_idx} = []; % Skip failed schedule
            end
        end
        
        % Remove empty entries and append to global results
        local_results = vertcat(local_results{:});
        results = [results; local_results]; %#ok<AGROW>
        
        % Increment backup counter
        backup_counter = backup_counter + numel(local_results);
        
        % Save backup if counter exceeds interval
        if mod(backup_counter, backup_interval) == 0
            backup_filename = fullfile(backup_folder, sprintf('laminate_results_backup_%d.mat', backup_counter));
            save(backup_filename, 'results');
            % fprintf('Backup saved: %s (Iterations: %d)\n', backup_filename, backup_counter);
        end
    end
end

% Final save of results
save('laminate_results.mat', 'results');
disp('Sweeps complete. Results saved to laminate_results.mat.');

% Clean up parallel pool
delete(gcp);
