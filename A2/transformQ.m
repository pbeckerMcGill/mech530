function [Q_transf] = transformQ(Q, theta)

% Convert theta from degrees to radians
theta_rad = deg2rad(theta);

% Extract components from Q matrix as per the provided layout
Qxx = Q(1, 1);
Qyy = Q(2, 2);
Qxy = Q(1, 2); % Assuming symmetry, Qxy = Qyx
Qss = Q(3, 3);

% Compute the trigonometric terms based on theta
cos_2theta = cos(2*theta_rad);
cos_4theta = cos(4*theta_rad);
sin_2theta = sin(2*theta_rad);
sin_4theta = sin(4*theta_rad);

% Calculate U1, U2, U3, U4, U5 using the given Q-matrix values
U1 = (1/8) * (3*Qxx + 3*Qyy + 2*Qxy + 4*Qss);
U2 = (1/2) * (Qxx - Qyy);
U3 = (1/8) * (Qxx + Qyy - 2*Qxy - 4*Qss);
U4 = (1/8) * (Qxx + Qyy + 6*Qxy - 4*Qss);
U5 = (1/8) * (Qxx + Qyy - 2*Qxy + 4*Qss);

% Rebuild the Q_transf matrix by calculating each element manually
Q11 = U1 * 1 + U2 * cos_2theta + U3 * cos_4theta;
Q22 = U1 * 1 - U2 * cos_2theta + U3 * cos_4theta;
Q12 = U4 * 1 - U3 * cos_4theta;
Q66 = U5 * 1 - U3 * cos_4theta;
Q16 = U2 * 0.5 * sin_2theta + U3 * sin_4theta;
Q26 = U2 * 0.5 * sin_2theta - U3 * sin_4theta;

% Rebuild the transformed Q matrix
Q_transf = [Q11, Q12, Q16;
            Q12, Q22, Q26;
            Q16, Q26, Q66];

end
