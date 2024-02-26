clear
clc
close all

% addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
load("data.txt")
data(:,3:4) = data(:,3:4)./1000; % convert to km

% constants 
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
r0_N = [757.700; 5222.607; 4851.500]; % [km]
v0_N = [2.21321; 4.67834; -5.3713]; % [km/s]
mu = 3.986004415e5; % [km^3/s^2] earth's gravitational parameter
J2 = 1.082626925638815e-3; % J2 perturbation
Cd = 2; % approximate coefficient of drag
R0s1_E = [-5127.5100; -3794.1600; 0.0]; % [km]
R0s2_E = [3860.9100; 3238.4900; 3898.0940]; % [km]
R0s3_E = [549.5050; -1380.8720; 6182.1970]; % [km]

X0 = [r0_N; v0_N; mu; J2; Cd; R0s1_E; R0s2_E; R0s3_E];
S0 = [X0; reshape(eye(n),[],1)];

% get time span
t = data(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,Sref] = ode45(@(t,S) dragJ2ODE(t,S,constants), t, S0, options);
% plotEarthOrbit(Sref(:,1:3)', ae, "Simulated Reference Trajectory")
Xref = Sref(:,1:n)';

dx0 = zeros(n,1);
P0 = diag([ones(6,1); 1e-61; 1e6; 1e6; 1e-16*ones(3,1); ones(6,1)]); 
noise_sd = [1e-5; 1e-6]; % km km/s
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];

max_iterations = 8;

%% CKF

[X_CKF, dx_CKF, P_CKF, y_CKF, alpha_CKF] = newCKF(data, R, X0, dx0, P0, 1, "both", constants);
titles = ["CKF Pre-Fit Residuals vs. Time", "CKF Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_CKF, alpha_CKF, noise_sd, titles)
RMSresidual_CKF = sqrt(1/length(t)*sum([alpha_CKF(1,:);alpha_CKF(2,:)].^2, 2));
RMSresidualpre_CKF = sqrt(1/length(t)*sum([y_CKF(1,:);y_CKF(2,:)].^2, 2));
sd_CKF = diag(P_CKF(:,:,end))'.^(1/2);
Xf_CKF = X_CKF(:,end)' + dx_CKF(:,end)';

trueX_CKF = X_CKF+dx_CKF;
deltaX_CKF = Xref -(X_CKF+dx_CKF);
plotDeltaX(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")
plotDeltaXsc(t, deltaX_CKF, "CKF $\delta x_{LS}$ vs Time")

[X_CKFiterated, dx_CKFiterated, P_CKFiterated, y_CKFiterated, alpha_CKFiterated] = newCKF(data, R, X0, dx0, P0, 20, "both", constants);
titles = ["Iterated CKF Pre-Fit Residuals vs. Time", "Iterated CKF Post-Fit Residuals vs. Time"];
plotResiduals(data(:,1), y_CKFiterated, alpha_CKFiterated, noise_sd, titles)
RMSresidual_iteratedCKF = sqrt(1/length(t)*sum([alpha_CKFiterated(1,:);alpha_CKFiterated(2,:)].^2, 2));
RMSresidualpre_iteratedCKF = sqrt(1/length(t)*sum([y_CKFiterated(1,:);y_CKFiterated(2,:)].^2, 2));
sd_CKFiterated = diag(P_CKFiterated(:,:,end))'.^(1/2);
Xf_CKFiterated = X_CKFiterated(:,end)' + dx_CKFiterated(:,end)';

CKFtrace = zeros(1,length(t));
CKFtrace_r = zeros(1,length(t));
CKFtrace_v = zeros(1,length(t));
CKFtraceIterated = zeros(1,length(t));
for i = 1:length(t)
    CKFtrace(i) = trace(P_CKF(1:6,1:6,i));
    CKFtrace_r(i) = trace(P_CKF(1:3,1:3,i));
    CKFtrace_v(i) = trace(P_CKF(4:6,4:6,i));
    CKFtraceIterated(i) = trace(P_CKFiterated(1:6,1:6,i));
end
figure
semilogy(t/60^2, CKFtrace)
title("Trace of the Covariance of Spacecraft State")
xlabel("time [hours]")
ylabel("trace(P)  [km^2 + (km/s)^2]")
figure
semilogy(t/60^2, CKFtraceIterated)
title("Trace of the Covariance of Spacecraft State (Iterated CKF)")
xlabel("time [hours]")
ylabel("trace(P)  [km^2 + (km/s)^2]")
% figure
% semilogy(t/60^2, CKFtrace_r)
% title("Trace of the Covariance of Position")
% xlabel("time [hours]")
% ylabel("trace(P) [(km/s)^2]")
% figure
% semilogy(t/60^2, CKFtrace_v)
% title("Trace of the Covariance of Velocity")
% xlabel("time [hours]")
% ylabel("trace(P) [km^2]")
%% Batch

[X_batch3, dx0_batch3, P_batch3, y_batch3, alpha_batch3, iterations_batch3, RMSresidual_batch3] = newbatch(data, R, X0, dx0, P0, 3, 'both', constants);
titles = ["Batch Pre-Fit Residuals vs. Time (Only 3 Iterations)", "Batch Post-Fit Residuals vs. Time (Only 3 Iterations)"];
% plotResiduals(data(:,1), y_batch3, alpha_batch3, noise_sd, titles)
% plotPostFitResiduals(data(:,1), alpha_batch3, noise_sd, titles)
plotResiduals(data(:,1), y_batch3, alpha_batch3, noise_sd, titles)
RMSresidualpre_3batch = sqrt(1/length(t)*sum([y_batch3(1,:);y_batch3(2,:)].^2, 2));

[X_batch, dx0_batch, P_batch, y_batch, alpha_batch, iterations_batch, RMSresidual_batch] = newbatch(data, R, X0, dx0, P0, max_iterations, 'both', constants);
titles = ["Batch Pre-Fit Residuals vs. Time", "Batch Post-Fit Residuals vs. Time"];
% plotResiduals(data(:,1), y_batch, alpha_batch, noise_sd, titles)
plotPostFitResiduals(data(:,1), alpha_batch, noise_sd, titles)
RMSresidualpre_batch = sqrt(1/length(t)*sum([y_batch(1,:);y_batch(2,:)].^2, 2));
sd_batch = diag(P_batch(:,:,end))'.^(1/2);
Xf_batch = X_batch(:,end)';

deltaX_batch = Xref - X_batch;
plotDeltaX(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
plotDeltaXsc(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")

%% Only one measurement type

% % CKF -------------------------------------------------------------------
% [X_CKF_range, dx_CKF_range, P_CKF_range, y_CKF_range, alpha_CKF_range] = newCKF(data, R, X0, dx0, P0, 1, "range", constants);
% titles = ["CKF Range Only Pre-Fit Residuals vs. Time", "CKF Range Only Post-Fit Residuals vs. Time"];
% plotRangeResiduals(data(:,1), y_CKF_range, alpha_CKF_range, noise_sd, titles)
% 
% [X_CKF_rangeRate, dx_CKF_rangeRate, P_CKF_rangeRate, y_CKF_rangeRate, alpha_CKF_rangeRate] = newCKF(data, R, X0, dx0, P0, 1, "rangeRate", constants);
% titles = ["CKF Range Rate Only Pre-Fit Residuals vs. Time", "CKF Range Rate Only Post-Fit Residuals vs. Time"];
% plotRangeRateResiduals(data(:,1), y_CKF_rangeRate, alpha_CKF_rangeRate, noise_sd, titles)
% 
% deltaX_CKF_range = Xref -(X_CKF_range+dx_CKF_range);
% plotDeltaX(t, deltaX_CKF, "Range Only CKF $\delta x_{LS}$ vs Time")
% deltaX_CKF_rangeRate = Xref -(X_CKF_range+dx_CKF_range);
% plotDeltaX(t, deltaX_CKF_rangeRate, "Range Rate Only CKF $\delta x_{LS}$ vs Time")
% 
% figure
% compareCovariance(P_CKF_range(1:3, 1:3, end), P_CKF_rangeRate(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["Range", "Range Rate"], "CKF $3\sigma$ Covariance Ellipses - Limited Measurements", 3.2e-5)

% Batch -------------------------------------------------------------------
[X_batch_range, dx_batch_range, P_batch_range, y_batch_range, alpha_batch_range, iterations_range, RMS_range] = newbatch(data, R, X0, dx0, P0, max_iterations, "range", constants);
titles = ["Batch Range Only Pre-Fit Residuals vs. Time", "Batch Range Only Post-Fit Residuals vs. Time"];
plotRangeResiduals(data(:,1), y_batch_range, alpha_batch_range, noise_sd, titles)

[X_batch_rangeRate, dx_batch_rangeRate, P_batch_rangeRate, y_batch_rangeRate, alpha_batch_rangeRate, iterations_range, RMS_rangeRate] = newbatch(data, R, X0, dx0, P0, max_iterations, "rangeRate", constants);
titles = ["Batch Range Rate Only Pre-Fit Residuals vs. Time", "Batch Range Rate Only Post-Fit Residuals vs. Time"];
plotRangeRateResiduals(data(:,1), y_batch_rangeRate, alpha_batch_rangeRate, noise_sd, titles)

figure
compareCovariance(P_batch_range(1:3, 1:3, end), P_batch_rangeRate(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["Range", "Range Rate"], "Batch $3\sigma$ Covariance Ellipses - Limited Measurements", 3.2e-5)
figure
compareCovariance(P_batch_range(4:6, 4:6, end), P_batch_rangeRate(4:6, 4:6, end), ["Vx [km/s]", "Vy [km/s]", "Vz [km/s]"], ["Range", "Range Rate"], "$3\sigma$ Covariance Ellipses of Velocity - Limited Measurement", 5.1e-8)
figure
compareCovariance(P_batch_range(1:3, 1:3, end), P_batch(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["Range", "Both"], "Batch $3\sigma$ Covariance Ellipses - Limited Measurements", 4.1e-5)

%% Changing fixed stations

P0 = diag([ones(6,1); 1e-61; 1e6; 1e6; ones(3,1); ones(6,1)]);

% [X_CKF_noFix, dx_CKF_noFix, P_CKF_noFix, y_CKF_noFix, alpha_CKF_noFix] = newCKF(data, R, X0, dx0, P0, 1, "both", constants);
% titles = ["No Fixed Station CKF Pre-Fit Residuals vs. Time", "No Fixed Station CKF Post-Fit Residuals vs. Time"];
% plotResiduals(data(:,1), y_CKF_noFix, alpha_CKF_noFix, noise_sd, titles)
% RMSresidual_CKF_noFix = sqrt(1/length(t)*sum([alpha_CKF_noFix(1,:);alpha_CKF_noFix(2,:)].^2, 2));
% deltaX_CKF_noFix = Xref -(X_CKF_noFix+dx_CKF_noFix);
% plotDeltaX(t, deltaX_CKF_noFix, "CKF $\delta x_{LS}$ vs Time")

[X_batch_noFix, dx0_batch_noFix, P_batch_noFix, y_batch_noFix, alpha_batch_noFix, iterations_batch_noFix, RMSresidual_batch_noFix] = newbatch(data, R, X0, dx0, P0, max_iterations, 'both', constants);
titles = ["No Fixed Station Batch Pre-Fit Residuals vs. Time", "No Fixed Station Batch Post-Fit Residuals vs. Time"];
plotPostFitResiduals(data(:,1), alpha_batch_noFix, noise_sd, titles)
% plotDeltaX(t, Xref-X_batch_noFix, "No Fixed Stations Batch $\delta x_{LS}$ vs Time")

% Fix Station 2 ------------------------------------------------------------
P0 = diag([ones(6,1); 1e-61; 1e6; 1e6; ones(3,1); 1e-16*ones(3,1); ones(3,1)]);

% [X_CKF_fix337, dx_CKF_fix337, P_CKF_fix337, y_CKF_fix337, alpha_CKF_fix337] = newCKF(data, R, X0, dx0, P0, 1, "both", constants);
% titles = ["337 Fixed CKF Pre-Fit Residuals vs. Time", "337 Fixed CKF Post-Fit Residuals vs. Time"];
% plotPostFitResiduals(data(:,1), alpha_CKF_fix337, noise_sd, titles)
% RMSresidual_CKF_fix337 = sqrt(1/length(t)*sum([alpha_CKF_fix337(1,:);alpha_CKF_fix337(2,:)].^2, 2));
% deltaX_CKF_fix337 = Xref -(X_CKF_fix337+dx_CKF_fix337);
% plotDeltaX(t, deltaX_CKF_fix337, "337 Fixed CKF $\delta x_{LS}$ vs Time")

[X_batch_fix337, dx0_batch_fix337, P_batch_fix337, y_batch_fix337, alpha_batch_fix337, iterations_batch_fix337, RMSresidual_batch_fix337] = newbatch(data, R, X0, dx0, P0, max_iterations, 'both', constants);
titles = ["337 Fixed Batch Pre-Fit Residuals vs. Time", "337 Fixed Batch Post-Fit Residuals vs. Time"];
plotPostFitResiduals(data(:,1), alpha_batch_fix337, noise_sd, titles)
% plotDeltaX(t, Xref-X_batch_fix337, "337 Fixed Batch $x_{LS}$ vs Time")

% Fix Station 3 ------------------------------------------------------------
P0 = diag([ones(6,1); 1e-61; 1e6; 1e6; ones(3,1); ones(3,1); 1e-16*ones(3,1)]);

% [X_CKF_fix394, dx_CKF_fix394, P_CKF_fix394, y_CKF_fix394, alpha_CKF_fix394] = newCKF(data, R, X0, dx0, P0, 1, "both", constants);
% titles = ["394 Fixed CKF Pre-Fit Residuals vs. Time", "394 Fixed CKF Post-Fit Residuals vs. Time"];
% plotPostFitResiduals(data(:,1), alpha_CKF_fix394, noise_sd, titles)
% RMSresidual_CKF_fix394 = sqrt(1/length(t)*sum([alpha_CKF_fix394(1,:);alpha_CKF_fix394(2,:)].^2, 2));
% deltaX_CKF_fix394 = Xref -(X_CKF_fix394+dx_CKF_fix394);
% plotDeltaX(t, deltaX_CKF_fix394, "394 Fixed CKF $\delta x_{LS}$ vs Time")

[X_batch_fix394, dx0_batch_fix394, P_batch_fix394, y_batch_fix394, alpha_batch_fix394, iterations_batch_fix394, RMSresidual_batch_fix394] = newbatch(data, R, X0, dx0, P0, max_iterations, 'both', constants);
titles = ["394 Fixed Batch Pre-Fit Residuals vs. Time", "394 Fixed Batch Post-Fit Residuals vs. Time"];
plotPostFitResiduals(data(:,1), alpha_batch_fix394, noise_sd, titles)
% plotDeltaX(t, Xref-X_batch_fix394, "394 Fixed Batch $\delta x_{LS}$ vs Time")

figure
compareCovariance(P_batch_noFix(1:3, 1:3, end), P_batch(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["None Fixed", "101 Fixed"], "", 6e-5)
compareCovariance(P_batch_fix337(1:3, 1:3, end), P_batch_fix394(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["None Fixed", "101 Fixed", "337 Fixed", "394 Fixed"], "$3\sigma$ Covariance Ellipses - Varying Fixed Station", 6e-5)

% deltaX_batch = Xref - X_batch_noFix;
% plotDeltaX(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
% % plotDeltaXsc(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
% deltaX_batch = Xref - X_batch_fix337;
% plotDeltaX(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
% % plotDeltaXsc(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
% deltaX_batch = Xref - X_batch_fix394;
% plotDeltaX(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
% % plotDeltaXsc(t, deltaX_batch, "Batch $\Delta x_{LS}$ vs Time")
%% Plotting
figure
compareCovariance(P_CKF(1:3, 1:3, end), P_batch(1:3, 1:3, end), ["X [km]", "Y [km]", "Z [km]"], ["CKF", "Batch"], "$3\sigma$ Covariance Ellipses of Position at $t_f$", 3.2e-5)
figure
compareCovariance(P_CKF(4:6, 4:6, end), P_batch(4:6, 4:6, end), ["Vx [km/s]", "Vy [km/s]", "Vz [km/s]"], ["CKF", "Batch"], "$3\sigma$ Covariance Ellipses of Velocity at $t_f$", 5.1e-8)
plotDeltaXsc(t, X_batch-trueX_CKF, "Batch - CKF vs Time")
plotDeltaX(t, X_batch-trueX_CKF, "Batch - CKF vs Time")