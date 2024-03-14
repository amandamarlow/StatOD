clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW1_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW3_6080')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\Project1')
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
Xideal = Sideal(:,1:n)';
t = t';

%% Setup

X0 = S0(1:n);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
% sigma_SNC = 10e-6;
sigma_SNC = 1e-6;
Qc = diag(sigma_SNC^2*ones(1,3));

%% Comparing with no process noise

[dx_smooth_noQ, state_smooth_noQ, P_smooth_noQ] = smoothedCKF(tideal, Y_simulated, R, zeros(3,3), X0, zeros(n,1), P0, constants);
smootherError_noQ = state_smooth_noQ - Xideal;
% plotErrorAndBounds_HW3(tspan, smootherError_noQ, P_smooth_noQ, "Smoothed CKF Error vs Time")

[X_batch, dx0_batch, P_batch, y_batch, alpha_batch, iterations_batch, RMSresidual_batch] = batch_HW4(Y_simulated, R, X0, zeros(n,1), P0, 1, constants);
batchError = X_batch-Xideal;

plotErrorAndBounds_HW4(tideal, smootherError_noQ, P_smooth_noQ, batchError, P_batch, "Smoothed CKF + Batch Error vs Time")

%% Comparing to batch with no J3

% j = 1;
% for i = 1:length(t)
%      tmeas = find(t(i) == YJ3_simulated(:,1), 1);
%      if ~isempty(tmeas)
%         tJ3_idx(j) = tmeas;
%         j = j+1;
%      end
% end
% 
% [dx_smooth_noQ, state_smooth_noQ, P_smooth_noQ] = smoothedCKF(YJ3_simulated(:,1), YJ3_simulated, R, zeros(3,3), X0, zeros(n,1), P0, constants);
% smootherError_noQ = state_smooth_noQ - Xtrue(:,tJ3_idx);
% % plotErrorAndBounds_HW3(tspan, smootherError_noQ, P_smooth_noQ, "Smoothed CKF Error vs Time")
% 
% [X_batch, dx0_batch, P_batch, y_batch, alpha_batch, iterations_batch, RMSresidual_batch] = batch_HW4(YJ3_simulated, R, X0, zeros(n,1), P0, 1, constants);
% batchError = X_batch-Xtrue(:,tJ3_idx);
% 
% plotErrorAndBounds_HW4(YJ3_simulated(:,1), smootherError_noQ, P_smooth_noQ, batchError, P_batch, "Smoothed CKF Error vs Time")

%% Comparing SNC to Smoothed with J3

[X_CKF_SNC, dx_CKF_SNC, P_CKF_SNC, y_CKF_SNC, alpha_CKF_SNC] = CKF_SNC(tspan, YJ3_simulated, R, Qc, X0, zeros(n,1), P0, constants);
error_SNC = (X_CKF_SNC + dx_CKF_SNC) - Xtrue;
plotErrorAndBounds_HW3(tspan, error_SNC, P_CKF_SNC, "CKF Error vs Time (SNC only)")

[dx_smooth, state_smooth, P_smooth] = smoothedCKF(tspan, YJ3_simulated, R, Qc, X0, zeros(n,1), P0, constants);
smootherError = state_smooth - Xtrue;
plotErrorAndBounds_HW3(tspan, smootherError, P_smooth, "Smoothed CKF Error vs Time")

RMSposition3d_SNC = rms(vecnorm(error_SNC(1:3,:),2,1), "omitnan");
RMSvelocity3d_SNC = rms(vecnorm(error_SNC(4:6,:),2,1), "omitnan");
RMSposition_SNC = rms(error_SNC(1:3,:), 2, "omitnan");
RMSvelocity_SNC = rms(error_SNC(4:6,:), 2, "omitnan");
RMSposition3d_smoother = rms(vecnorm(smootherError(1:3,:),2,1), "omitnan");
RMSvelocity3d_smoother = rms(vecnorm(smootherError(4:6,:),2,1), "omitnan");
RMSposition_smoother = rms(smootherError(1:3,:),2, "omitnan");
RMSvelocity_smoother = rms(smootherError(4:6,:),2, "omitnan");
