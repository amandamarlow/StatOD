clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW3_6080\')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
load("YJ3_simulated.mat")
load("initialConds.mat")
load("constants.mat")
load("XtrueJ3.mat")

% % given orbit elements w/ respect to ECI
% a = 10000; %km
% e = 0.001; 
% i = 40*pi/180; %rad
% OMEGA = 80*pi/180; %rad
% omega = 40*pi/180; %rad
% nu0 = 0;
% 
% % constants
% mu = 398600.4415; % [km^3/s^2] earth's gravitational parameter
% constants.mu = mu;
% ae = 6378.0; % [km] mean equitorial radius of earth
% constants.ae = ae; % [km] mean equitorial radius of earth
% J2 = 0.0010826269; % J2 perturbation
% constants.J2 = J2;
J3 = -0.0000025323;
constants.J3 = J3;
% C = [mu; J2; J3];
% 
% constants.omegaE = 7.2921158553e-5; % [rad/s]
% constants.theta0 = 122*pi/180; % [rad]
% c = 299792.458; % [km/s] speed of light
% 
% %% Simulate True Trajectory
% 
% n = 6; % length of state vector
% 
% % convert orbit elements to initial state in ECI frame
% rp = a*(1-e);
% r0_p = [rp; 0; 0];
% vp = sqrt((1+e)/(1-e)*mu/a);
% v0_p = [0; vp; 0];
% PN = M3(omega)*M1(i)*M3(OMEGA);
% r0_ECI = PN'*r0_p;
% v0_ECI = PN'*v0_p;
% S0 = [r0_ECI; v0_ECI; C; reshape(eye(9),[],1)];
% 
% % % get time span
% % T = 2*pi*sqrt(a^3/mu);
% % tspan = 0:10:15*T;
% % load("Y_simulated.mat")
% % tspanIdeal = Y_simulated(:,1);
% % load("YJ3_simulated.mat")
% % % tspan = YJ3_simulated(:,1);
% 
% % simulate of reference trajectory
% options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% % [tideal,Sideal] = ode45(@(t,S) orbitODE(t,S,constants), tspanIdeal, S0, options);
% [t,S] = ode45(@(t,S) orbitJ3ODE(t,S,constants), YJ3_simulated(:,1)', S0, options);
% XtrueJ3 = S(:,1:n)';
% t = t';


%% Implement Filters with ideal sigma
n = 6;
% X0 = S0(1:n);
X0 = XtrueJ3(:,1);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
% P0 = blkdiag(0.001*eye(3), (1e-6)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
% dx0 = [0.5, 0.5, 0.5, 0.5e-3, 0.5e-3, 0.5e-3]';
dx0 = [ones(1,3)*5, 5*ones(1,3)*10^-3]';

sigma_SNC = 1e-6;
Qc = diag(sigma_SNC^2*ones(1,3));

warmStart = 1;
% [X_EKF, P_EKF, preFit_EKF, postFit_EKF] = EKF_SNC(tspan, YJ3_simulated, R, Qc, X0, P0, warmStart, constants);
% error_EKF = X_EKF - XtrueJ3;
% [X_EKF_dx0, P_EKF_dx0, preFit_EKF_dx0, postFit_EKF_dx0] = EKF_SNC(tspan, YJ3_simulated, R, Qc, X0+dx0, P0, warmStart, constants);
% error_EKF_dx0 = X_EKF_dx0 - XtrueJ3;

alpha = 1;
beta = 2;
% [x_UKF, P_UKF, preFit_UKF] = UKF(tspan, YJ3_simulated, R, zeros(3), X0, P0, alpha, beta, constants);
% error_UKF = x_UKF - XtrueJ3;
% 
% [x_UKF_Q, P_UKF_Q, preFit_UKF_Q] = UKF(tspan, YJ3_simulated, R, Qc, X0, P0, alpha, beta, constants);
% error_UKF_Q = x_UKF_Q - XtrueJ3;

% [x_UKF_dx0, P_UKF_dx0, preFit_UKF_dx0] = UKF(tspan, YJ3_simulated, R, Qc, X0+dx0, P0, alpha, beta, constants);
% error_UKF_dx0 = x_UKF_dx0 - XtrueJ3;

% alpha = 10^-4;
% [x_UKF_alpha2, P_UKF_alpha2, preFit_UKF_alpha2] = UKF(tspan, YJ3_simulated, R, Qc, X0, P0, alpha, beta, constants);
% error_UKF_alpha2 = x_UKF_alpha2 - XtrueJ3;

alpha = 1;
[x_UKF_J3, P_UKF_J3, preFit_UKF_J3] = UKF_J3(tspan, YJ3_simulated, R, zeros(3), X0, P0, alpha, beta, constants);
error_UKF_J3 = x_UKF_J3 - XtrueJ3;

%% RMS

% RMSposition_EKF_SNC = rms(vecnorm(error_EKF_SNC(1:3,:),2,1), "omitnan");
% RMSvelocity_EKF_SNC = rms(vecnorm(error_EKF_SNC(4:6,:),2,1), "omitnan");
% RMSresiduals_EKF_SNC = rms(alpha_EKF_SNC, 2, "omitnan");

%% Plot 

% % titles = ["EKF Pre-Fit Residuals vs. Time", "EKF-SNC Post-Fit Residuals vs. Time"];
% % plotResiduals_HW3(tspan, preFit_EKF, postFit_EKF, noise_sd, titles)
% plotErrorAndBounds_HW3(tspan, error_EKF, P_EKF, "EKF-SNC Error vs Time")

% titles = ["EKF Pre-Fit Residuals vs. Time Perturbed", "EKF-SNC Post-Fit Residuals vs. Time Perturbed"];
% plotResiduals_HW3(tspan, preFit_EKF, postFit_EKF, noise_sd, titles)
% plotErrorAndBounds_HW3(tspan, error_EKF_dx0, P_EKF_dx0, "EKF-SNC Error vs Time (Perturbed)")

% % titles = ["UKF Pre-Fit Residuals vs. Time", "UKF Post-Fit Residuals vs. Time"];
% plotErrorAndBounds_HW3(tspan, error_UKF, P_UKF, "UKF Error vs Time")
% % plotErrorAndBounds_HW3(YJ3_simulated(:,1)', error_UKF, P_UKF, "UKF Error vs Time")

% % % titles = ["UKF Pre-Fit Residuals vs. Time with Process Noise", "UKF Post-Fit Residuals vs. Time with Process Noise"];
% plotErrorAndBounds_HW3(tspan, error_UKF_Q, P_UKF_Q, "UKF Error vs Time with Process Noise")

% plotErrorAndBounds_HW3(tspan, error_UKF_dx0, P_UKF_dx0, "UKF Error vs Time (Perturbed)")
% plotErrorAndBounds_HW3(tspan, error_UKF_alpha2, P_UKF_alpha2, "UKF Error vs Time ($\alpha=1e-4$)")
plotErrorAndBounds_HW3(tspan, error_UKF_J3, P_UKF_J3, "UKF Error vs Time with J3")