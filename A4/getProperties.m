function [E_x, E_y, E_s, nu_x, nu_y, m, X_t, X_c, Y_t, Y_c, S_c, h_o, rho] = getProperties(pathDB, compositeName)
%GETPROPERTIESFROMDB Summary of this function goes here
%   Detailed explanation goes here

jsonText = fileread(pathDB);
jsonData = jsondecode(jsonText);

E_x = jsonData.(compositeName).Modulus.Ex;
E_y = jsonData.(compositeName).Modulus.Ey;
E_s = jsonData.(compositeName).Modulus.Es;
nu_x = jsonData.(compositeName).Modulus.vx;
nu_y.value = nu_x.value * (E_y.value / E_x.value);
nu_y.unit = "none";

m.value = (1 - nu_x.value*nu_y.value)^(-1);
m.unit = "none";

X_t = jsonData.(compositeName).Strength.Xt;
X_c = jsonData.(compositeName).Strength.Xc;
Y_t = jsonData.(compositeName).Strength.Yt;
Y_c = jsonData.(compositeName).Strength.Yc;
S_c = jsonData.(compositeName).Strength.Sc;

h_o = jsonData.(compositeName).Other.ho;
rho = jsonData.(compositeName).Other.rho;
end

