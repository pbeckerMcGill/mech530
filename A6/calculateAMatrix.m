function [A] = calculateAMatrix(layup_angles, h_i,  Q)

% Q's are converted to Pa
Qxx = Q(1, 1) * 1e6;
Qyy = Q(2, 2) * 1e6;
Qxy = Q(1, 2) * 1e6;
Qss = Q(3, 3) * 1e6;

% U's are converted to Pa
U1 = (1/8) * (3*Qxx + 3*Qyy + 2*Qxy + 4*Qss);
U2 = (1/2) * (Qxx - Qyy);
U3 = (1/8) * (Qxx + Qyy - 2*Qxy - 4*Qss);
U4 = (1/8) * (Qxx + Qyy + 6*Qxy - 4*Qss);
U5 = (1/8) * (Qxx + Qyy - 2*Qxy + 4*Qss);

% Get the number of plies
num_plies = length(layup_angles);

% Initialize the V_star values
V_star_1 = 0;
V_star_2 = 0;
V_star_3 = 0;
V_star_4 = 0;

% Loop through each ply angle
for i = 1:num_plies
    theta = layup_angles(i);
    
    theta_2 = deg2rad(2 * theta);
    theta_4 = deg2rad(4 * theta);
    
    V_star_1 = V_star_1 + cos(theta_2) * h_i;
    V_star_2 = V_star_2 + cos(theta_4) * h_i;
    V_star_3 = V_star_3 + sin(theta_2) * h_i;
    V_star_4 = V_star_4 + sin(theta_4) * h_i;
end

% Divide each V_star by the number of plies
V1 = 2 * V_star_1;
V2 = 2 * V_star_2;
V3 = 2 * V_star_3;
V4 = 2 * V_star_4;

h = 2 * num_plies * h_i;

A11 = h * U1 + U2 * V1 + U3 * V2;
A22 = h * U1 - U2 * V1 + U3 * V2;
A12 = h * U4 - U3 * V2;
A66 = h * U5 - U3 * V2;
A16 = 0.5 * U2 * V3 + U3 * V4;
A26 = 0.5 * U2 * V3 - U3 * V4;

A = [A11, A12, A16;
     A12, A22, A26;
     A16, A26, A66];

end
