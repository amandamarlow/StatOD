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
% J3 = -0.0000025323;
J3 = 0;
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
[t1,S1] = ode45(@(t,S) orbitODE(t,S,constants), tspan, S0, options);

% STM
flatSTM = S1(:,10:end);
STM = reshape(flatSTM',9,9,[]);

%% Simulate Measurements

% station locations
latlon_s1 = [-35.39833; 148.981944]*pi/180; % [rad]
latlon_s2 = [40.427222; 355.749444]*pi/180; % [rad]
latlon_s3 = [35.247164; 243.205]*pi/180; % [rad]
r_s1_E = latlon2ECEF(ae, latlon_s1(1), latlon_s1(2));
r_s2_E = latlon2ECEF(ae, latlon_s2(1), latlon_s2(2));
r_s3_E = latlon2ECEF(ae, latlon_s3(1), latlon_s3(2));

[range_observations, rangeRate_observations, elevations_byStation, elevations_all] = simMeas(t1, S1', [r_s1_E,r_s2_E,r_s3_E], constants);

% Generate Data 
Y_ideal = [range_observations(:,1:2), range_observations(:,3), rangeRate_observations(:,3)]; % [t(s), station#, range measurement, range rate measurement]
% with noise
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
range_noise = mvnrnd(0,(noise_sd(1))^2, length(range_observations));
rangeRate_noise = mvnrnd(0,(noise_sd(2))^2, length(rangeRate_observations));
Y_onlynoise = [range_observations(:,1:2), range_observations(:,3) + range_noise, rangeRate_observations(:,3) + rangeRate_noise]; % [t(s), station#, range measurement, range rate measurement]
idx1 = find(Y_onlynoise(:,2)==1);
idx2 = find(Y_onlynoise(:,2)==2);
idx3 = find(Y_onlynoise(:,2)==3);

%% Implement Filters

% MAKE SURE STATE INCLUDES CONSTANTS (ITS ALREADY IN THE ORBITODE/INTEGRATION
% FUNCTIONS)

% % Organize measurements for filters
% %   Y = Measurement Matrix -> [time after epoch [s], station number, measurement type
% numR = length(rangeMeasurements);
% numRR = length(rangeRateMeasurements);
% Y = cell(numR+numRR,4);
% R = zeros(numR+numRR,1);
% q = 1;
% for i = 1:max([numR, numRR])
%    if i <= numR
%        Y(q,1) = {range_observations(i,1)};
%        Y(q,2) = {range_observations(i,2)};
%        Y(q,3) = {"range"};
%        Y(q,4) = {range_observations(i,3)};
%        R(q) = (1e-3)^2;
%        q = q+1;
%    elseif i <= numRR
%        Y(q,1) = {rangeRate_observations(i,1)};
%        Y(q,2) = {rangeRate_observations(i,2)};
%        Y(q,3) = {"rangeRate"};
%        Y(q,4) = {rangeRate_observations(i,3)};
%        R(q) = (1e-6)^2;
%        q = q+1;
%    end
% end

n = 6; % length of state vector
m = length(Y_onlynoise);
X0 = S0(1:n);
dx0 = zeros(n,1);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
R = zeros(2,2,m);
R(1,1,:) = (noise_sd(2))^2;
R(2,2,:) = (noise_sd(2))^2;

% CKF - Classic/Linearized Kalman Filter
[t, X_CKF_ideal, dx_CKF_ideal, P_CKF_ideal, y_CKF_ideal, alpha_CKF_ideal] = CKF(Y_ideal, R, X0, dx0, P0, constants);
% error_CKF_ideal = (X_CKF_ideal + dx_CKF_ideal) - S1(:,1:3)';
[t, X_CKF_onlynoise, dx_CKF_onlynoise, P_CKF_onlynoise, y_CKF_onlynoise, alpha_CKF_onlynoise] = CKF(Y_onlynoise, R, X0, dx0, P0, constants);
% error_CKF_onlynoise = (X_CKF_onlynoise + dx_CKF_onlynoise) - S1(:,1:3)';
%% Plotting

% plot simulated measurements
figure
subplot(2,1,1)
sgtitle("Simulated Measurements")
scatter(Y_onlynoise(idx1, 1), Y_onlynoise(idx1, 3), '.')
hold on
scatter(Y_onlynoise(idx2, 1), Y_onlynoise(idx2, 3), '.')
scatter(Y_onlynoise(idx3, 1), Y_onlynoise(idx3, 3), '.')
legend('Station 1', 'Station 2', 'Station 3', 'Location', 'northeast')
xlabel('Time [s]')
ylabel("Range ($\rho$) [km]", 'Interpreter', 'latex')
grid on
subplot(2,1,2)
scatter(Y_onlynoise(idx1, 1), Y_onlynoise(idx1, 4), '.')
hold on
scatter(Y_onlynoise(idx2, 1), Y_onlynoise(idx2, 4), '.')
scatter(Y_onlynoise(idx3, 1), Y_onlynoise(idx3, 4), '.')
legend('Station 1', 'Station 2', 'Station 3', 'Location', 'northeast')
xlabel('Time [s]')
ylabel("Range Rate ($\dot{\rho}$) [km/s]", 'Interpreter', 'latex')
grid on

% Plot Residuals
titles = ["Pre-Fit Residuals vs. Time (no perturbation + no noise)", "Post-Fit Residuals vs. Time (no perturbation + no noise)"];
plotResiduals(t, y_CKF_ideal, alpha_CKF_ideal, noise_sd, titles)
titles = ["Pre-Fit Residuals vs. Time (no perturbation + noise)", "Post-Fit Residuals vs. Time (no perturbation + noise)"];
plotResiduals(t, y_CKF_onlynoise, alpha_CKF_onlynoise, noise_sd, titles)

% Plot Error vs time and 3 sigma bounds
% plotErrorAndBounds(error_CKF_ideal, P_CKF_ideal)
% plotErrorAndBounds(error_CKF_onlynoise, P_CKF_onlynoise)