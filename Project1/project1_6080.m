clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
% addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW1_6080')
% addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
% load("data.txt")
% data(:,3:4) = data(:,3:4)./1000; % convert to km
data = load("project.txt");
data(:,3:4) = data(:,3:4)/1000; % convert to km

% constants 
% ae = 6378136.3; % [m] mean equitorial radius of earth
% constants.ae = ae; % [m] mean equitorial radius of earth
% area = 3; % m^2
% constants.area = area;
% rho0 = 3.614e-13; % kg/m^3
% constants.rho0 = rho0;
% r0 = 700000.0 + ae; % m
% constants.r0 = r0;
% H0 = 88667.0; %m
% constants.H0 = H0;
% m = 970; % kg
% constants.m = m;
% constants.omegaE = 7.2921158553e-5; % [rad/s]
% constants.theta0 = 0; % [rad]
ae = 6378136.3e-3; % [km] mean equitorial radius of earth
constants.ae = ae; % [km] mean equitorial radius of earth
area = 3e-6; % km^2
constants.area = area;
rho0 = 3.614e-4; % kg/km^3
constants.rho0 = rho0;
r0 = 700000.0e-3 + ae; % km
constants.r0 = r0;
H0 = 88667.0e-3; % km
constants.H0 = H0;
m = 970; % kg
constants.m = m;
constants.omegaE = 7.2921158553e-5; % [rad/s]
constants.theta0 = 0; % [rad]


%% Simulate True Trajectory

n = 18; % length of state vector
% Initial Conditions
% r0_N = [757700; 5222607; 4851500]; % [m]
r0_N = [757.700; 5222.607; 4851.500]; % [km]
% v0_N = [2213.21; 4678.34; -5371.3]; % [m/s]
v0_N = [2.21321; 4.67834; -5.3713]; % [km/s]
% mu = 3.986004415e14; % [m^3/s^2] earth's gravitational parameter
mu = 3.986004415e5; % [km^3/s^2] earth's gravitational parameter
J2 = 1.082626925638815e-3; % J2 perturbation
Cd = 2; % approximate coefficient of drag
% R0s1_E = [-5127510.0; -3794160.0; 0.0]; % [m]
% R0s2_E = [3860910.0; 3238490.0; 3898094.0]; % [m]
% R0s3_E = [549505.0; -1380872.0; 6182197.0]; % [m]
R0s1_E = [-5127.5100; -3794.1600; 0.0]; % [km]
R0s2_E = [3860.9100; 3238.4900; 3898.0940]; % [km]
R0s3_E = [549.5050; -1380.8720; 6182.1970]; % [km]

X0 = [r0_N; v0_N; mu; J2; Cd; R0s1_E; R0s2_E; R0s3_E];
S0 = [X0; reshape(eye(n),[],1)];

% get time span
% t = data(1,1):20:data(end,1);
t = data(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,Sref] = ode45(@(t,S) dragJ2ODE(t,S,constants), t, S0, options);
% plotEarthOrbit(Sref(:,1:3)', ae, "Simulated Reference Trajectory")
Xref = Sref(:,1:n)';

dx0 = zeros(n,1);
% P0 = diag([1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e20, 1e6, 1e6, 1e-10, 1e-10, 1e-10, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6]);
P0 = diag([ones(6,1); 1e-61; 1e6; 1e6; 1e-16*ones(3,1); ones(6,1)]); 
noise_sd = [1e-5; 1e-6]; % km km/s
% R = zeros(2,2,length(data(:,1)));
% R(1,1,:) = (noise_sd(1))^2;
% R(2,2,:) = (noise_sd(2))^2;
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];

%% CKF

[X_CKF, dx_CKF, P_CKF, y_CKF, alpha_CKF] = newCKF(data, R, X0, dx0, P0, 1, constants);
titles = ["CKF Pre-Fit Residuals vs. Time", "CKF Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_CKF, alpha_CKF, noise_sd, titles)

RMSresidual_CKF = sqrt(1/length(t)*sum([alpha_CKF(1,:);alpha_CKF(2,:)].^2, 2));

deltaX_CKF = Xref - (X_CKF+dx_CKF);
% plotDeltaX(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")
plotDeltaXsc(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")

CKFtrace = zeros(1,length(t));
for i = 1:length(t)
    CKFtrace(i) = trace(P_CKF(1:6,1:6,i));
end
figure
semilogy(t/60^2, CKFtrace)
title("Trace of the Covariance of Spacecraft State")
xlabel("time [hours]")
ylabel("trace(P)  [km^2] & [(km/s)^2]")

[X_CKFiterated, dx_CKFiterated, P_CKFiterated, y_CKFiterated, alpha_CKFiterated] = newCKF(data, R, X0, dx0, P0, 20, constants);
titles = ["Iterated CKF Pre-Fit Residuals vs. Time", "Iterated CKF Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_CKFiterated, alpha_CKFiterated, noise_sd, titles)
deltaX_CKF_iterated = Xref - (X_CKF+dx_CKF);
% plotDeltaX(t, deltaX_CKF_iterated, "CKF $\delta x_{LS}$ vs Time")
plotDeltaXsc(t, deltaX_CKF_iterated, "CKF $\delta x_{LS}$ vs Time")
%% Batch

[X_batch, dx0_batch, P_batch, y_batch, alpha_batch, iterations_batch, RMSresidual_batch] = newbatch(data, R, X0, dx0, P0, constants);
titles = ["Batch Pre-Fit Residuals vs. Time", "Batch Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_batch, alpha_batch, noise_sd, titles)

deltaX_batch = Xref - X_batch;
% plotDeltaX(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
plotDeltaXsc(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
