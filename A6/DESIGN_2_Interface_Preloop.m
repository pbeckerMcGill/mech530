clc;

% Inputs
material = "flaxpreg";
[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", material);
chosen_angle = 55;
schedule = [chosen_angle, -chosen_angle, chosen_angle, -chosen_angle, chosen_angle, -chosen_angle, chosen_angle, -chosen_angle, chosen_angle, -chosen_angle, chosen_angle, -chosen_angle];
schedule_full = [schedule, flip(schedule)];
symmetrical = true;

z_c = 0; % m

% Compute S's
S_xx = 1 / E_x.value;
S_xy = -nu_y.value / E_y.value;
S_yx = -nu_x.value / E_x.value;
S_yy = 1 / E_y.value;
S_ss = 1 / E_s.value;

% Compute Q's
Q_xx = m.value * E_x.value;
Q_yy = m.value * E_y.value;
Q_yx = m.value * nu_x.value * E_y.value;
Q_xy = m.value * nu_y.value * E_x.value;
Q_ss = E_s.value;

% Build matrices
S = [S_xx S_xy 0; S_yx S_yy 0; 0 0 S_ss];
Q = [Q_xx Q_xy 0; Q_yx Q_yy 0; 0 0 Q_ss];

fprintf('CHOSEN MATERIAL: %s\n\n', material);

fprintf("\n================= MATRICES =================\n");

D_matrix = calculateDMatrix(schedule, h_o.value, Q, z_c);
d_matrix = inv(D_matrix);

fprintf('\n[D] Matrix (in %s):\n', "Nm");
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(1, 1), D_matrix(1, 2), D_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(2, 1), D_matrix(2, 2), D_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', D_matrix(3, 1), D_matrix(3, 2), D_matrix(3, 3));

fprintf('\n[d] Matrix (in %s):\n', "(Nm)^-1");
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(1, 1), d_matrix(1, 2), d_matrix(1, 3));
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(2, 1), d_matrix(2, 2), d_matrix(2, 3));
fprintf('%15.3e %15.3e %15.3e\n', d_matrix(3, 1), d_matrix(3, 2), d_matrix(3, 3));