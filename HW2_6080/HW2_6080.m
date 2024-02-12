clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
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

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,S] = ode45(@(t,S) orbitODE(t,S,constants), tspan, S0, options);


%% Simulate Measurements

% station locations
latlon_s1 = [-35.39833; 148.981944]*pi/180; % [rad]
latlon_s2 = [40.427222; 355.749444]*pi/180; % [rad]
latlon_s3 = [35.247164; 243.205]*pi/180; % [rad]
r_s1_E = latlon2ECEF(ae, latlon_s1(1), latlon_s1(2));
r_s2_E = latlon2ECEF(ae, latlon_s2(1), latlon_s2(2));
r_s3_E = latlon2ECEF(ae, latlon_s3(1), latlon_s3(2));

[range_observations, rangeRate_observations, elevations_byStation, elevations_all] = simMeas(t, S', [r_s1_E,r_s2_E,r_s3_E], constants);

% Generate Data 
Y_ideal = [range_observations(:,1:2), range_observations(:,3), rangeRate_observations(:,3)]; % [t(s), station#, range measurement, range rate measurement]
% with noise
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
% range_noise = mvnrnd(0,(noise_sd(1))^2, length(range_observations));
% rangeRate_noise = mvnrnd(0,(noise_sd(2))^2, length(rangeRate_observations));
% Y_simulated = [range_observations(:,1:2), range_observations(:,3) + range_noise, rangeRate_observations(:,3) + rangeRate_noise]; % [t(s), station#, range measurement, range rate measurement]
load("Y_simulated.mat")
idx1 = find(Y_simulated(:,2)==1);
idx2 = find(Y_simulated(:,2)==2);
idx3 = find(Y_simulated(:,2)==3);

%% Implement Filters

n = 6; % length of state vector
m = length(Y_simulated);
X0 = S0(1:n);
dx0_ideal = zeros(n,1);
dx0 = [0.9; 0.9; 0.9; 0.8e-3; 0.8e-3; 0.8e-3];
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
R = zeros(2,2,m);
R(1,1,:) = (noise_sd(1))^2;
R(2,2,:) = (noise_sd(2))^2;


%% CKF - Classic/Linearized Kalman Filter

% [X_CKF_ideal, dx_CKF_ideal, P_CKF_ideal, y_CKF_ideal, alpha_CKF_ideal] = CKF(t, Y_ideal, R, X0, dx0_ideal, P0, constants);
% error_CKF_ideal = (X_CKF_ideal + dx_CKF_ideal) - S(:,1:n)';
% [X_CKF_onlynoise, dx_CKF_onlynoise, P_CKF_onlynoise, y_CKF_onlynoise, alpha_CKF_onlynoise] = CKF(t, Y_simulated, R, X0, dx0_ideal, P0, constants);
% error_CKF_onlynoise = (X_CKF_onlynoise + dx_CKF_onlynoise) - S(:,1:n)';
[X_CKF_dx0noise, dx_CKF_dx0noise, P_CKF_dx0noise, y_CKF_dx0noise, alpha_CKF_dx0noise] = CKF(t, Y_simulated, R, X0, dx0, P0, constants);
error_CKF_dx0noise = (X_CKF_dx0noise + dx_CKF_dx0noise) - S(:,1:n)';

% get RMS
dataTime = t>Y_simulated(1,1);
RMSerror_CKF_component = sqrt(1/length(t)*sum(error_CKF_dx0noise.^2, 2));
RMSerror_CKF_component_afterPass = sqrt(1/sum(dataTime)*sum(error_CKF_dx0noise(:,dataTime).^2, 2));
RMSerror_CKF_3D = sqrt(1/length(t)*sum([norm(error_CKF_dx0noise(1:3,:));norm(error_CKF_dx0noise(4:6,:))].^2, 2));
RMSerror_CKF_3D_afterPass = sqrt(1/sum(dataTime)*sum([norm(error_CKF_dx0noise(1:3,dataTime));norm(error_CKF_dx0noise(4:6,dataTime))].^2, 2));
RMSresidual_CKF = sqrt(1/length(t)*sum(alpha_CKF_dx0noise.^2, 2));
RMSresidual_CKF_afterPass = sqrt(1/sum(dataTime)*sum(alpha_CKF_dx0noise(:,dataTime).^2, 2));

% Plot Residuals
% titles = ["CKF Pre-Fit Residuals vs. Time (Ideal)", "CKF Post-Fit Residuals vs. Time (Ideal)"];
% plotResiduals(t, y_CKF_ideal, alpha_CKF_ideal, noise_sd, titles)
% titles = ["CKF Pre-Fit Residuals vs. Time (noise)", "CKF Post-Fit Residuals vs. Time (noise)"];
% plotResiduals(t, y_CKF_onlynoise, alpha_CKF_onlynoise, noise_sd, titles)
titles = ["CKF Pre-Fit Residuals vs. Time ($\delta x_0$ + noise)", "CKF Post-Fit Residuals vs. Time ($\delta x_0$ + noise)"];
plotResiduals(t, y_CKF_dx0noise, alpha_CKF_dx0noise, noise_sd, titles)

% Plot Error vs time and 3 sigma bounds
% plotErrorAndBounds(t, error_CKF_ideal, P_CKF_ideal, "CKF Error vs Time (Ideal)")
% plotErrorAndBounds(t, error_CKF_onlynoise, P_CKF_onlynoise, "CKF Error vs Time (no perturbation + noise)")
plotErrorAndBounds(t, error_CKF_dx0noise, P_CKF_dx0noise, "CKF Error vs Time ($\delta x_0$ + noise)")


%% EKF - Extended Kalman Filter

% % [X_EKF_ideal, P_EKF_ideal, y_EKF_ideal, alpha_EKF_ideal] = EKF(t, Y_ideal, R, X0, P0, constants);
% % error_EKF_ideal = X_EKF_ideal - S(:,1:n)';
% % [X_EKF_onlynoise, P_EKF_onlynoise, y_EKF_onlynoise, alpha_EKF_onlynoise] = EKF(t, Y_simulated, R, X0, P0, constants);
% % error_EKF_onlynoise = X_EKF_onlynoise - S(:,1:n)';
[X_EKF_dx0noise, P_EKF_dx0noise, y_EKF_dx0noise, alpha_EKF_dx0noise] = EKF(t, Y_simulated, R, X0+dx0, P0, constants);
error_EKF_dx0noise = X_EKF_dx0noise - S(:,1:n)';

% get RMS
dataTime = t>Y_simulated(1,1);
RMSerror_EKF_component = sqrt(1/length(t)*sum(error_EKF_dx0noise.^2, 2));
RMSerror_EKF_component_afterPass = sqrt(1/sum(dataTime)*sum(error_EKF_dx0noise(:,dataTime).^2, 2));
RMSerror_EKF_3D = sqrt(1/length(t)*sum([norm(error_EKF_dx0noise(1:3,:));norm(error_EKF_dx0noise(4:6,:))].^2, 2));
RMSerror_EKF_3D_afterPass = sqrt(1/sum(dataTime)*sum([norm(error_EKF_dx0noise(1:3,dataTime));norm(error_EKF_dx0noise(4:6,dataTime))].^2, 2));
RMSresidual_EKF = sqrt(1/length(t)*sum(alpha_EKF_dx0noise.^2, 2));
RMSresidual_EKF_afterPass = sqrt(1/sum(dataTime)*sum(alpha_EKF_dx0noise(:,dataTime).^2, 2));

% % Plot Residuals
% titles = ["EKF Pre-Fit Residuals vs. Time (Ideal)", "EKF Post-Fit Residuals vs. Time (Ideal)"];
% plotResiduals(t, y_EKF_ideal, alpha_EKF_ideal, noise_sd, titles)
% titles = ["EKF Pre-Fit Residuals vs. Time (noise)", "EKF Post-Fit Residuals vs. Time (noise)"];
% plotResiduals(t, y_EKF_onlynoise, alpha_EKF_onlynoise, noise_sd, titles)
titles = ["EKF Pre-Fit Residuals vs. Time ($\delta x_0 + noise$)", "EKF Post-Fit Residuals vs. Time ($\delta x_0$ + noise)"];
plotResiduals(t, y_EKF_dx0noise, alpha_EKF_dx0noise, noise_sd, titles)
% 
% % Plot Error vs time and 3 sigma bounds
% plotErrorAndBounds(t, error_EKF_ideal, P_EKF_ideal, "EKF Error vs Time (Ideal)")
% plotErrorAndBounds(t, error_EKF_onlynoise, P_EKF_onlynoise, "EKF Error vs Time (noise)")
plotErrorAndBounds(t, error_EKF_dx0noise, P_EKF_dx0noise, "EKF Error vs Time ($\delta x_0 + noise$)")


%% Batch Filter

% [X_batch_ideal, dx_batch_ideal, P_batch_ideal, y_batch_ideal, alpha_batch_ideal] = batch(t, Y_ideal, R, X0, dx0_ideal, P0, constants);
% error_batch_ideal = X_batch_ideal - S(:,1:n)';
% [X_batch_onlynoise, dx_batch_onlynoise, P_batch_onlynoise, y_batch_onlynoise, alpha_batch_onlynoise] = batch(t, Y_simulated, R, X0, dx0_ideal, P0, constants);
% error_batch_onlynoise = X_batch_onlynoise - S(:,1:n)';
% [X_batch_dx0noise, dx_batch_dx0noise, P_batch_dx0noise, y_batch_dx0noise, alpha_batch_dx0noise] = batch(t, Y_simulated, R, X0, dx0_ideal, P0, constants);
% error_batch_dx0noise = X_batch_dx0noise - S(:,1:n)';
% 
% % Plot Residuals
% titles = ["batch Pre-Fit Residuals vs. Time (Ideal)", "batch Post-Fit Residuals vs. Time (Ideal)"];
% plotResiduals(t, y_batch_ideal, alpha_batch_ideal, noise_sd, titles)
% titles = ["batch Pre-Fit Residuals vs. Time (noise)", "batch Post-Fit Residuals vs. Time (noise)"];
% plotResiduals(t, y_batch_onlynoise, alpha_batch_onlynoise, noise_sd, titles)
% titles = ["batch Pre-Fit Residuals vs. Time ($\delta x_0 + noise$)", "batch Post-Fit Residuals vs. Time ($\delta x_0$ + noise)"];
% plotResiduals(t, y_batch_dx0noise, alpha_batch_dx0noise, noise_sd, titles)
% 
% % Plot Error vs time and 3 sigma bounds
% plotErrorAndBounds(t, error_batch_ideal, P_batch_ideal, "batch Error vs Time (Ideal)")
% plotErrorAndBounds(t, error_batch_onlynoise, P_batch_onlynoise, "batch Error vs Time (noise)")
% plotErrorAndBounds(t, error_batch_dx0noise, P_batch_dx0noise, "batch Error vs Time ($\delta x_0 + noise$)")


%% Other Plotting

% % plot simulated measurements
% figure
% subplot(2,1,1)
% sgtitle("Simulated Measurements")
% scatter(Y_simulated(idx1, 1), Y_simulated(idx1, 3), '.')
% hold on
% scatter(Y_simulated(idx2, 1), Y_simulated(idx2, 3), '.')
% scatter(Y_simulated(idx3, 1), Y_simulated(idx3, 3), '.')
% legend('Station 1', 'Station 2', 'Station 3', 'Location', 'northeast')
% xlabel('Time [s]')
% ylabel("Range ($\rho$) [km]", 'Interpreter', 'latex')
% grid on
% subplot(2,1,2)
% scatter(Y_simulated(idx1, 1), Y_simulated(idx1, 4), '.')
% hold on
% scatter(Y_simulated(idx2, 1), Y_simulated(idx2, 4), '.')
% scatter(Y_simulated(idx3, 1), Y_simulated(idx3, 4), '.')
% legend('Station 1', 'Station 2', 'Station 3', 'Location', 'northeast')
% xlabel('Time [s]')
% ylabel("Range Rate ($\dot{\rho}$) [km/s]", 'Interpreter', 'latex')
% grid on
