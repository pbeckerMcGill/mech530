function [strain_vector_transf] = transformStrain(strain_vector, theta)
% Convert theta from degrees to radians
theta_rad = deg2rad(theta);

% Calculate m and n
m = cos(theta_rad);
n = sin(theta_rad);

% Define the transformation matrix using m and n
T = [m^2, n^2, m*n;
    n^2, m^2, -m*n;
    -2*m*n, 2*m*n, m^2 - n^2];

% Perform the transformation by multiplying the transformation matrix T
% with the input strain vector
strain_vector_transf = T * strain_vector;
end