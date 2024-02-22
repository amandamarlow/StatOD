function [Sdot] = dragJ2ODE(t,S,constants)
%ORBITODE Summary of this function goes here
%   S is the 6x1 state vector [position; velocity]

ae = constants.ae;
omegaE = constants.omegaE;
omegaE_N = [0;0;omegaE];
theta0 = constants.theta0;
area = constants.area;
rho0 = constants.rho0;
r0 = constants.r0;
H0 = constants.H0;
m = constants.m;

n = 18; % length of state vector
n_const = n-6;

r_N = S(1:3);
r = norm(r_N);
v_N = S(4:6);

mu = S(7);
J2 = S(8);
Cd = S(9);
Rs1_E = S(10:12);
Rs2_E = S(13:15);
Rs3_E = S(16:18);

Phi = reshape(S(n+1:end),[n,n]);

% DCM for ECI and ECEF frames
alpha = t*omegaE + theta0;
EN = Euler3(alpha);
NE = EN';

%% All of this section is in ECEF coordinates
r_E = EN*r_N;
x = r_E(1);
y = r_E(2);
z = r_E(3);

% Dynamics
a_J2_E = -3*mu*(ae^2)*J2/2/(r^7) * [x*(r^2 - 5*z^2); y*(r^2 - 5*z^2); z*(r^2 + 2*(x^2 + y^2) - 3*z^2)];

% Jacobians
aJ2_E_partial_R_E = [
    -((3 * J2 * ae^2 * (-7*x^2*(x^2+y^2-4*z^2) + r^2*(3*x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2))), (3 * J2 * ae^2 * x*y * (-2*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * x*z * (8*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2));
    (3 * J2 * ae^2 * x*y * (-2*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2)), -((3 * J2 * ae^2 * (-7*y^2*(x^2+y^2-4*z^2) + r^2*(x^2+3*y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2))), (3 * J2 * ae^2 * y*z * (8*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2));
    (3 * J2 * ae^2 * x*z * (-6*r^2 + 21*(x^2+y^2) - 14*z^2) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * y*z * (-6*r^2 + 21*(x^2+y^2) - 14*z^2) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * (21*(x^2+y^2)*z^2 - 14*z^4 - 3*r^2*(x^2+y^2-2*z^2)) * mu) / (2 * r^8 * sqrt(r^2))
];
aJ2_E_partial_C = [
    -((3 * J2 * ae^2 * x * (x^2+y^2-4*z^2)) / (2 * (x^2+y^2+z^2)^(7/2))), -((3 * ae^2 * x * (x^2+y^2-4*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2))), zeros(1,n_const-2);
    -((3 * J2 * ae^2 * y * (x^2+y^2-4*z^2)) / (2 * (x^2+y^2+z^2)^(7/2))), -((3 * ae^2 * y * (x^2+y^2-4*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2))), zeros(1,n_const-2);
    (3 * J2 * ae^2 * z * (-3*(x^2+y^2)+2*z^2)) / (2 * (x^2+y^2+z^2)^(7/2)), (3 * ae^2 * z * (-3*(x^2+y^2)+2*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2)), zeros(1,n_const-2)
];

%% This section in in ECI coordinates

x = r_N(1);
y = r_N(2);
z = r_N(3);

% drag
vr_N = v_N - cross(omegaE_N,r_N); % atmospheric relative velocity
rho = rho0 * exp(-(r-r0)/H0); % exponential density model
a_D_N = -1/2*rho*Cd*area/m*norm(vr_N)*vr_N;
%gravity
a_mu_N = -mu*r_N/(r^3);
a_J2_N = NE*a_J2_E;
%total acceleration
a_N = a_mu_N + a_J2_N + a_D_N;

% Jacobians
amu_Partial_R_N = [
    -((sqrt(r^2) * (r^2 - 3*x^2) * mu) / r^6), (3*x*y*mu) / (r^2)^(5/2), (3*x*z*mu) / (r^2)^(5/2);
    (3*x*y*mu) / (r^2)^(5/2), -((sqrt(r^2) * (r^2 - 3*y^2) * mu) / r^6), (3*y*z*mu) / (r^2)^(5/2);
    (3*x*z*mu) / (r^2)^(5/2), (3*y*z*mu) / (r^2)^(5/2), -((sqrt(r^2) * (r^2 - 3*z^2) * mu) / r^6)
];
% amu_Partial_C = [
%     -(x/r^3), 0, 0;
%     -(y/r^3), 0, 0;
%     -(z/r^3), 0, 0
% ];
amu_Partial_C = [-r_N/r^3, zeros(3, n_const-1)];
aD_N_partial_R_N = Cd*area*rho/2/m * (norm(vr_N)*(vr_N*r_N'/H0/r + tilde(omegaE_N)) + vr_N*vr_N'/norm(vr_N)*tilde(omegaE_N));
aD_N_partial_V_N = -rho*Cd*area/2/m*(vr_N*vr_N'/norm(vr_N) + norm(vr_N)*eye(3));
aD_N_partial_C = [zeros(3,2), -rho*area/2/m*norm(vr_N)*vr_N, zeros(3,9)];

aJ2_N_partial_R_N = NE*aJ2_E_partial_R_E*EN;
aJ2_N_partial_C = NE*aJ2_E_partial_C;


% STM Integration
a_partial_C = amu_Partial_C + aJ2_N_partial_C + aD_N_partial_C;
a_partial_R_N = amu_Partial_R_N + aJ2_N_partial_R_N + aD_N_partial_R_N;
A = [
    zeros(3), eye(3), zeros(3, n_const);
    a_partial_R_N, aD_N_partial_V_N, a_partial_C;
    zeros(n_const, 3), zeros(n_const, 3), zeros(n_const);
];


PhiDot = A*Phi;

Sdot = [v_N; a_N; zeros(n-6,1); reshape(PhiDot,[],1)];
end

