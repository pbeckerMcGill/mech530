name_material = "graphite_epoxy_1";
[E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties("material_database.json", name_material);
material_properties = [E_x.value, E_y.value, E_s.value, nu_x.value, nu_y.value, m.value, X_t.value, X_c.value, Y_t.value, Y_c.value, S_c.value, h_o.value, rho.value];
z_c = 0.005; % m
schedule_base = [0, 0, 20, -20, 0, 90];
symmetrical = true;

P_1 = 1000; % [N]
L = 0.5; % [m]
b = 0.1; % [m]
M_1 = -(P_1 * L) / (4 * b);
N_1 = (0.5 * P_1) / b;

M_vector = [M_1; 0; 0];
N_vector = [N_1; 0; 0];