function [epsilon_i_aboveCore_innerSide, epsilon_i_aboveCore_outerSide, epsilon_i_belowCore_innerSide, epsilon_i_belowCore_outerSide] = calculateStrain(layup_angles, h_i, z_c, epsilon_o_vector, k_vector)

% Computes epsilon_i from epsilon_i = epsilon_i_o + z * k_i

% Get the number of plies
num_plies = length(layup_angles);

z_0 = z_c;
z_i1 = z_0;

z_aboveCore_innerSide = zeros(1, num_plies);
z_aboveCore_outerSide = zeros(1, num_plies);

% Loop through each ply angle, above the core
for i = 1:num_plies
    % Calculate z_i and z_(i-1)
    z_i = z_i1 + h_i; % Increment z by h_i for each layer

    % Append elements to list
    z_aboveCore_innerSide(i) = z_i1;
    z_aboveCore_outerSide(i) = z_i;
    
    % Update z_i1 for the next iteration
    z_i1 = z_i;
end


z_belowCore_innerSide = -z_aboveCore_innerSide;
z_belowCore_outerSide = -z_aboveCore_outerSide;

epsilon_i_o = num_plies * epsilon_o_vector;

epsilon_i_aboveCore_innerSide = epsilon_i_o + sum(z_aboveCore_innerSide) * k_vector;
epsilon_i_aboveCore_outerSide = epsilon_i_o + sum(z_aboveCore_outerSide) * k_vector;
epsilon_i_belowCore_innerSide = epsilon_i_o + sum(z_belowCore_innerSide) * k_vector;
epsilon_i_belowCore_outerSide = epsilon_i_o + sum(z_belowCore_outerSide) * k_vector;

end