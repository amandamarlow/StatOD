clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW1_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW3_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\Project1')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
load("Y_simulated.mat")
load("YJ3_simulated.mat")
load("constants.mat")
load("S0.mat")

%% Householder check
% eq 5.6.8
A = [1, -2, -1; 2, -1, 1; 1, 1, 2];
[Q, TA] = qr(A);

% residuals won't really give more information

% simulate reference trajectory
tspan = Y_simulated(1,1):10:Y_simulated(end,1);
n = 6;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,S] = ode45(@(t,S) orbitODE(t,S,constants), tspan, S0, options);
[tJ3,SJ3] = ode45(@(t,S) orbitJ3ODE(t,S,constants), tspan, S0, options);
Xtrue = S(:,1:n)';
XtrueJ3 = SJ3(:,1:n)';
t = t';

% Initial Condition
X0 = S0(1:n);
dx0 = [0.5; 0.5; 0.5; 0.5e-3; 0.5e-3; 0.5e-3];
% dx0 = zeros(n,1);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
meas_cov = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];

%% Question 1

[X_basic, P_basic, Xref_basic, dx_basic, ytilde_basic, e_basic] = SRIF_basic(t, Y_simulated, X0+dx0, zeros(n,1), P0, meas_cov, constants);
error_basic = X_basic-Xtrue;

[X_CKF, dx_CKF, P_CKF, y_CKF, alpha_CKF] = CKF_SNC(t, Y_simulated, meas_cov, zeros(3), X0+dx0, zeros(n,1), P0, constants);
error_CKF = (X_CKF + dx_CKF) - Xtrue;

plotErrorAndBounds(tspan, error_basic, P_basic, "Basic SRIF Error vs Time")
plotErrorAndBounds_HW5(tspan, error_basic, P_basic, error_CKF, P_CKF, "CKF + Basic SRIF Error vs Time")

%% Question 2
sigma_SNC = 1e-6;
Qc = diag(sigma_SNC^2*ones(1,3));

[X_procNoise, P_procNoise, Xref_procNoise, dx_procNoise, ytilde_procNoise, e_procNoise] = SRIF_procNoise(t, YJ3_simulated, X0, dx0, P0, meas_cov, Qc, constants);

[X_SNC, dx_SNC, P_SNC, y_SNC, alpha_SNC] = CKF_SNC(t, YJ3_simulated, meas_cov, Qc, X0+dx0, zeros(n,1), P0, constants);
error_SNC = (X_SNC + dx_SNC) - XtrueJ3;
% plotErrorAndBounds(tspan, error_SNC, P_SNC, "LKF SNC Error vs Time")
