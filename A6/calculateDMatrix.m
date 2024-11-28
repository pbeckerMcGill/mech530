function [D] = calculateDMatrix(layup_angles, h_i, Q, z_c)

layup_angles_bottom_up = flip(layup_angles);

% Q's are converted to Pa
Qxx = Q(1, 1) * 1e6;
Qyy = Q(2, 2) * 1e6;
Qxy = Q(1, 2) * 1e6;
Qss = Q(3, 3) * 1e6;

U1 = (1/8) * (3*Qxx + 3*Qyy + 2*Qxy + 4*Qss);
U2 = (1/2) * (Qxx - Qyy);
U3 = (1/8) * (Qxx + Qyy - 2*Qxy - 4*Qss);
U4 = (1/8) * (Qxx + Qyy + 6*Qxy - 4*Qss);
U5 = (1/8) * (Qxx + Qyy - 2*Qxy + 4*Qss);

% Get the number of plies
num_plies = length(layup_angles_bottom_up);

% Initialize the V_star values
V_star_1 = 0;
V_star_2 = 0;
V_star_3 = 0;
V_star_4 = 0;

% Initialize z_0
z_0 = z_c;
z_i1 = z_0;

% Display initial z_0
% fprintf('Initial z_0: %f\n', z_0);

% Loop through each ply angle
for i = 1:num_plies
    theta = layup_angles_bottom_up(i);
    
    % Calculate z_i and z_(i-1)
    z_i = z_i1 + h_i; % Increment z by h_i for each layer
    
    % Display z_i1, z_i, and the associated layup angle
    % fprintf('Iteration %d:\n', i);
    % fprintf('  z_i1: %f\n', z_i1);
    % fprintf('  z_i: %f\n', z_i);
    % fprintf('  Layup angle (theta): %f degrees\n', theta);
    
    % Angle calculations
    theta_2 = deg2rad(2 * theta);
    theta_4 = deg2rad(4 * theta);
    
    % Update V_star using the difference z_i^3 - z_(i-1)^3
    V_star_1 = V_star_1 + cos(theta_2) * (z_i^3 - z_i1^3);
    V_star_2 = V_star_2 + cos(theta_4) * (z_i^3 - z_i1^3);
    V_star_3 = V_star_3 + sin(theta_2) * (z_i^3 - z_i1^3);
    V_star_4 = V_star_4 + sin(theta_4) * (z_i^3 - z_i1^3);
    
    % Update z_i1 for the next iteration
    z_i1 = z_i;
end

h = 2 * (num_plies * h_i + z_c);
z_star_c = (2 * z_c) / h;
h_star = ((1 - z_star_c^3) * h^3) / 12;

V_star_1 = 8 / ((h^3) * (1 - z_star_c^3)) * V_star_1;
V_star_2 = 8 / ((h^3) * (1 - z_star_c^3)) * V_star_2;
V_star_3 = 8 / ((h^3) * (1 - z_star_c^3)) * V_star_3;
V_star_4 = 8 / ((h^3) * (1 - z_star_c^3)) * V_star_4;

% Calculate D-matrix components
D_star_11 = U1 + U2 * V_star_1 + U3 * V_star_2; % BAD
D_star_22 = U1 - U2 * V_star_1 + U3 * V_star_2; % BAD
D_star_12 = U4 - U3 * V_star_2; % GOOD
D_star_66 = U5 - U3 * V_star_2; % GOOD
D_star_16 = 0.5 * U2 * V_star_3 + U3 * V_star_4; % SIGN
D_star_26 = 0.5 * U2 * V_star_3 - U3 * V_star_4; % SIGN

D_star = [D_star_11, D_star_12, D_star_16;
          D_star_12, D_star_22, D_star_26;
          D_star_16, D_star_26, D_star_66];

D = D_star * h_star;

end
