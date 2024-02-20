clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
% addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW1_6080')
% addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
load("data.txt")

% constants 
ae = 6378.0e3; % [m] mean equitorial radius of earth
constants.ae = ae; % [m] mean equitorial radius of earth
area = 3; % m^2
constants.area = area;
rho0 = 3.614e-13; % kg/m^3
constants.rho0 = rho0;
r0 = 700000 + ae; % m
constants.r0 = r0;
H0 = 88667; %m
constants.H0 = H0;
m = 970; % kg
constants.m = m;
constants.omegaE = 7.2921158553e-5; % [rad/s]
constants.theta0 = 0; % [rad]
% c = 299792.458e3; % [m/s] speed of light

%% Simulate True Trajectory

n = 18; % length of state vector
% Initial Conditions
r0_N = [757700; 5222607; 4851500]; % [m]
v0_N = [2213.21; 4678.34; -5371.3]; % [m/s]
mu = 398600.4415e9; % [m^3/s^2] earth's gravitational parameter
J2 = 0.0010826269; % J2 perturbation
Cd = 2; % approximate coefficient of drag
R0s1_E = [-5127510.0; -3794160.0; 0.0]; % [m]
R0s2_E = [3860910.0; 3238490.0; 3898094.0]; % [m]
R0s3_E = [549505.0; -1380872.0; 6182197.0]; % [m]

X0 = [r0_N; v0_N; mu; J2; Cd; R0s1_E; R0s2_E; R0s3_E];
S0 = [X0; reshape(eye(n),[],1)];

% get time span
tspan = data(1,1):20:data(end,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,S] = ode45(@(t,S) dragJ2ODE(t,S,constants), tspan, S0, options);
plotEarthOrbit(S(:,1:3)', ae, "Simulated Reference Trajectory")

dx0 = zeros(n,1);
P0 = diag([1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, 1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6]);
% P0 = diag([1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e10, 1e6, 1e6, 1e-1, 1e-1, 1e-1, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6]);
noise_sd = [0.01; 0.001]; % m m/s
R = zeros(2,2,m);
R(1,1,:) = (noise_sd(1))^2;
R(2,2,:) = (noise_sd(2))^2;

[X_CKF, dx_CKF, P_CKF, y_CKF, alpha_CKF] = CKF(t, data, R, X0, dx0, P0, constants);
titles = ["CKF Pre-Fit Residuals vs. Time", "CKF Post-Fit Residuals vs. Time"];
plotResiduals(t, y_CKF, alpha_CKF, noise_sd, titles)

[X_batch, dx0_batch, P_batch, y_batch, alpha_batch, iterations_batch, RMSresidual_batch] = batch(t, data, R, X0, dx0, P0, constants);
titles = ["Batch Pre-Fit Residuals vs. Time", "Batch Post-Fit Residuals vs. Time"];
plotResiduals(t, y_batch, alpha_batch, noise_sd, titles)
