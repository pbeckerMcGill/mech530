function [S_transf] = transformS(S, theta)

% Convert theta from degrees to radians
theta_rad = deg2rad(theta);

% Extract components from the S matrix as per the provided layout
Sxx = S(1, 1);
Syy = S(2, 2);
Sxy = S(1, 2); % Assuming symmetry, Sxy = Syx
Sss = S(3, 3);

% Compute the trigonometric terms based on theta
cos_2theta = cos(2*theta_rad);
cos_4theta = cos(4*theta_rad);
sin_2theta = sin(2*theta_rad);
sin_4theta = sin(4*theta_rad);

% Calculate U1, U2, U3, U4, U5 using the given S-matrix values
U1 = (1/8) * (3*Sxx + 3*Syy + 2*Sxy + Sss);
U2 = (1/2) * (Sxx - Syy);
U3 = (1/8) * (Sxx + Syy - 2*Sxy - Sss);
U4 = (1/8) * (Sxx + Syy + 6*Sxy - Sss);
U5 = (1/2) * (Sxx + Syy - 2*Sxy + Sss);

% Rebuild the S_transf matrix by calculating each element manually
S11 = U1 * 1 + U2 * cos_2theta + U3 * cos_4theta;
S22 = U1 * 1 - U2 * cos_2theta + U3 * cos_4theta;
S12 = U4 * 1 - U3 * cos_4theta;
S66 = U5 * 1 - U3 * 4 * cos_4theta;
S16 = U2 * sin_2theta + U3 * 2 * sin_4theta;
S26 = U2 * sin_2theta - U3 * 2 * sin_4theta;

% Rebuild the transformed S matrix
S_transf = [S11, S12, S16;
            S12, S22, S26;
            S16, S26, S66];

end
