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


%% Implement Filters with ideal sigma
n = 6;
% X0 = S0(1:n);
X0 = XtrueJ3(:,1);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
% P0 = blkdiag(0.001*eye(3), (1e-6)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
dx0 = [0.5, 0.5, 0.5, 0.5e-3, 0.5e-3, 0.5e-3]';

sigma_SNC = 1e-6;
Qc = diag(sigma_SNC^2*ones(1,3));

warmStart = 1;
[X_EKF, P_EKF, preFit_EKF, postFit_EKF] = EKF_SNC(tspan, YJ3_simulated, R, Qc, X0, P0, warmStart, constants);
error_EKF = X_EKF - XtrueJ3;
[X_EKF_dx0, P_EKF_dx0, preFit_EKF_dx0, postFit_EKF_dx0] = EKF_SNC(tspan, YJ3_simulated, R, Qc, X0+dx0, P0, warmStart, constants);
error_EKF_dx0 = X_EKF_dx0 - XtrueJ3;

alpha = 1;
beta = 2;
[x_UKF, P_UKF, preFit_UKF] = UKF(tspan, YJ3_simulated, R, zeros(3), X0, P0, alpha, beta, constants);
error_UKF = x_UKF - XtrueJ3;

[x_UKF_Q, P_UKF_Q, preFit_UKF_Q] = UKF(tspan, YJ3_simulated, R, Qc, X0, P0, alpha, beta, constants);
error_UKF_Q = x_UKF_Q - XtrueJ3;

[x_UKF_dx0, P_UKF_dx0, preFit_UKF_dx0] = UKF(tspan, YJ3_simulated, R, Qc, X0+dx0, P0, alpha, beta, constants);
error_UKF_dx0 = x_UKF_dx0 - XtrueJ3;

alpha = 10^-4;
[x_UKF_alpha2, P_UKF_alpha2, preFit_UKF_alpha2] = UKF(tspan, YJ3_simulated, R, Qc, X0, P0, alpha, beta, constants);
error_UKF_alpha2 = x_UKF_alpha2 - XtrueJ3;

%% RMS

% RMSposition_EKF_SNC = rms(vecnorm(error_EKF_SNC(1:3,:),2,1), "omitnan");
% RMSvelocity_EKF_SNC = rms(vecnorm(error_EKF_SNC(4:6,:),2,1), "omitnan");
% RMSresiduals_EKF_SNC = rms(alpha_EKF_SNC, 2, "omitnan");

%% Plot 

% titles = ["EKF Pre-Fit Residuals vs. Time", "EKF-SNC Post-Fit Residuals vs. Time"];
        % plotResiduals_HW3(tspan, preFit_EKF, postFit_EKF, noise_sd, titles)
plotErrorAndBounds_HW3(tspan, error_EKF, P_EKF, "EKF-SNC Error vs Time")

% titles = ["EKF Pre-Fit Residuals vs. Time Perturbed", "EKF-SNC Post-Fit Residuals vs. Time Perturbed"];
% plotResiduals_HW3(tspan, preFit_EKF, postFit_EKF, noise_sd, titles)
plotErrorAndBounds_HW3(tspan, error_EKF_dx0, P_EKF_dx0, "EKF-SNC Error vs Time")

% titles = ["UKF Pre-Fit Residuals vs. Time", "UKF Post-Fit Residuals vs. Time"];
plotErrorAndBounds_HW3(tspan, error_UKF, P_UKF, "UKF Error vs Time")

% titles = ["UKF Pre-Fit Residuals vs. Time with Process Noise", "UKF Post-Fit Residuals vs. Time with Process Noise"];
plotErrorAndBounds_HW3(tspan, error_UKF_Q, P_UKF_Q, "UKF Error vs Time with Process Noise")

plotErrorAndBounds_HW3(tspan, error_UKF_dx0, P_UKF_dx0, "UKF Error vs Time perturbed")
plotErrorAndBounds_HW3(tspan, error_UKF_alpha2, P_UKF_alpha2, "UKF Error vs Time ($\alpha=1e-4$)")