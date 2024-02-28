clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW1_6080')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")

% given orbit elements w/ respect to ECI
a = 10000; %km
e = 0.001; 
i = 40*pi/180; %rad
OMEGA = 80*pi/180; %rad
omega = 40*pi/180; %rad
nu0 = 0;

% constants
mu = 398600.4415; % [km^3/s^2] earth's gravitational parameter
constants.mu = mu;
ae = 6378.0; % [km] mean equitorial radius of earth
constants.ae = ae; % [km] mean equitorial radius of earth
J2 = 0.0010826269; % J2 perturbation
constants.J2 = J2;
J3 = -0.0000025323;
C = [mu; J2; J3];

constants.omegaE = 7.2921158553e-5; % [rad/s]
constants.theta0 = 122*pi/180; % [rad]
c = 299792.458; % [km/s] speed of light

%% Simulate True Trajectory

n = 6; % length of state vector

% convert orbit elements to initial state in ECI frame
rp = a*(1-e);
r0_p = [rp; 0; 0];
vp = sqrt((1+e)/(1-e)*mu/a);
v0_p = [0; vp; 0];
PN = M3(omega)*M1(i)*M3(OMEGA);
r0_ECI = PN'*r0_p;
v0_ECI = PN'*v0_p;
S0 = [r0_ECI; v0_ECI; C; reshape(eye(9),[],1)];

% get time span
T = 2*pi*sqrt(a^3/mu);
tspan = 0:10:15*T;
load("Y_simulated.mat")
tspanIdeal = Y_simulated(:,1);
load("YJ3_simulated.mat")
% tspan = YJ3_simulated(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tideal,Sideal] = ode45(@(t,S) orbitODE(t,S,constants), tspanIdeal, S0, options);
[t,S] = ode45(@(t,S) orbitJ3ODE(t,S,constants), tspan, S0, options);
Xtrue = S(:,1:n)';
t = t';

%% Implement Filters
n = 6;
X0 = S0(1:n);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
warmStart = 1;

% sigma = 10e-15;
sigma = 1e-6;
Qc = diag(sigma^2*ones(1,3));

% [X_CKF_SNC, dx_CKF_SNC, P_CKF_SNC, y_CKF_SNC, alpha_CKF_SNC] = CKF_SNC(t, YJ3_simulated, R, Qc, X0, zeros(n,1), P0, constants);
% [X_EKF_SNC, P_EKF_SNC, y_EKF_SNC, alpha_EKF_SNC] = EKF_SNC(t, YJ3_simulated, R, Qc, X0, P0, warmStart, constants);

n = 9;
tau = T/30;
sigma_DMC = 1e-6;
Qc = diag(sigma_DMC^2*ones(1,3));
P0_aDMC = diag(ones(3,1) * 3*(sigma_DMC)^2*tau/2);
a_DMC = zeros(3,1);
X0_DMC = [X0; a_DMC];
P0_DMC = blkdiag(P0, P0_aDMC);

% [X_CKF_DMC, dx_CKF_DMC, P_CKF_DMC, y_CKF_DMC, alpha_CKF_DMC] = CKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, zeros(n,1), P0_DMC, tau, constants);
[X_EKF_DMC, P_EKF_DMC, y_EKF_DMC, alpha_EKF_DMC] = EKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, P0_DMC, tau, warmStart, constants);

%% Plot - SNC
% titles = ["CKF-SNC Pre-Fit Residuals vs. Time", "CKF-SNC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_CKF_SNC, alpha_CKF_SNC, noise_sd, titles)
% error_CKF_SNC = (X_CKF_SNC + dx_CKF_SNC) - Xtrue;
% plotErrorAndBounds_HW3(t, error_CKF_SNC, P_CKF_SNC, "CKF-SNC Error vs Time")
% 
% titles = ["EKF Pre-Fit Residuals vs. Time (J3)", "EKF-SNC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_EKF_SNC, alpha_EKF_SNC, noise_sd, titles)
% error_EKF_SNC = X_EKF_SNC - Xtrue;
% plotErrorAndBounds_HW3(t, error_EKF_SNC, P_EKF_SNC, "EKF-SNC Error vs Time")

%% Plot - DMC

% titles = ["CKF-DMC Pre-Fit Residuals vs. Time", "CKF-DMC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_CKF_DMC, alpha_CKF_DMC, noise_sd, titles)
% error_CKF_DMC = (X_CKF_DMC(1:6,:)+ dx_CKF_DMC(1:6,:)) - Xtrue;
% plotErrorAndBounds_HW3(t, error_CKF_DMC, P_CKF_DMC, "CKF-DMC Error vs Time")

titles = ["EKF-DMC Pre-Fit Residuals vs. Time", "EKF-DMC Post-Fit Residuals vs. Time"];
plotResiduals_HW3(t, y_EKF_DMC, alpha_EKF_DMC, noise_sd, titles)
error_EKF_DMC = X_EKF_DMC(1:6,:) - Xtrue;
plotErrorAndBounds_HW3(t, error_EKF_DMC, P_EKF_DMC, "EKF-DMC Error vs Time")