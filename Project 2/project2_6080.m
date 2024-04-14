clear
% clc
close all

addpath('C:\Users\marlo\OneDrive - UCB-O365\Documents\6080 - Statistical Orbit Determination\Project 2\Provided Files')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')

% readtable("Project2a_Obs_modified.txt")
data_2a = readmatrix("Project2a_Obs.txt");
data_2a(:,end) = [];
% readtable('Project2 Prob2 truth traj 50days.txt')
truthTraj = readmatrix("Project2_Prob2_truth_traj_50days.txt");
truthTraj(:,end) = [];

% setup synamics
constants.t0 = 2456296.25; % JD
constants.muS = 132712440017.987; % km^3/s^2 gravitiational parameter of sun
% constants.muE = 398600.4415; % [km^3/s^2] earth's gravitational parameter
srp_flux = 1357; % W/m^2 at 1 AU
c = 299792458; % m/s
constants.AU = 149597870.7; % km
constants.p_srAU = srp_flux/c*1000; % kN
constants.a2m = 0.01e-6; % km^2/kg


% measurement setup
constants.theta0 = 0;
constants.omegaE = 7.29211585275553e-5; % rad/s
constants.ae = 6378.1363; % km
DSS_34 = [-35.398333*pi/180, 148.981944*pi/180, 0.691750]; % [lat rad; lon rad; altitude km]
DSS_65 = [40.427222*pi/180, 355.749444*pi/180, 0.834539]; % [lat rad; lon rad; altitude km] NEGATIVE ONLY IN PART 2
DSS_13 = [35.247164*pi/180, 243.205*pi/180, 1.07114904]; % [lat rad; lon rad; altitude km]
constants.DSS_latlon = [DSS_34(1:2); DSS_65(1:2); DSS_13(1:2)];
constants.DSS_alt = [DSS_34(3); DSS_65(3); DSS_13(3)];

rngNoise = 5e-3; % km
rngRtNoise = 0.5e-6; % km
meas_cov = [rngNoise^2, 0; 0, rngRtNoise^2];

%% Estimate the State with a Known Target and Models

X0_true = truthTraj(1,2:8)'; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
n = length(X0_true);
sigX = 100; % km
sigV = 0.1; % km/s
sigCr = 0.1;
P0_2a = diag([sigX^2*ones(1,3), sigV^2*ones(1,3), sigCr^2]);
S0_traj = [X0_true; reshape(eye(n),[],1)];

% compare to provided truth trajectory
tspan_traj = truthTraj(:,1);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_traj,S_traj] = ode45(@(t,S) project2ODE(t,S,constants), tspan_traj, S0_traj, options);
r_traj = S_traj(:,1:3)';

figure
plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:))
hold on
% plotBody([0;0;0], constants.ae*1000)
plot3(truthTraj(:,2), truthTraj(:,3), truthTraj(:,4), '--')
axis equal
grid on
hold off

trajError = truthTraj-[t_traj, S_traj];
final_error = trajError(end, :);
fprintf("final state errors:")
disp(final_error(1:8))
fprintf("final STM errors: \n")
disp(reshape(final_error(9:end), n,n))

figure
subplot(2,1,1)
semilogy(abs(trajError))
subplot(2,1,2)
semilogy(vecnorm(abs(trajError),2,2))

% simulate whole trajectory
as = 696000; % km
tspan_2a = data_2a(:,1);
X0_2a = [-274096790.0; -92859240.0; -40199490.0; 32.67; -8.94; -3.88; 1.2]; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
n = length(X0_2a);
% S0_2a = [X0_2a; reshape(eye(n),[],1)];
S0_true = [X0_true; reshape(eye(n),[],1)];
[t_2a,S_2a] = ode45(@(t,S) project2ODE(t,S,constants), tspan_2a, S0_true, options);
X_2a_true = S_2a(:,1:n)';

colorS = "#EDB120";
colorE = "#77AC30";

r_S2E_N = zeros(3,length(tspan_2a));
for i = 1:length(tspan_2a)
    [r_S2E_N(:,i), ~, ~] = Ephem(constants.t0+tspan_2a(i)/86400,3, 'EME2000');
end
r_E2S_N = -r_S2E_N;

figure
plot3(X_2a_true(1,:), X_2a_true(2,:), X_2a_true(3,:))
hold on
plotBody([0;0;0], constants.ae*200, colorE)
plot3(r_E2S_N(1,:), r_E2S_N(2,:), r_E2S_N(3,:), 'Color', colorS)
plotBody(-r_S2E_N(:,end), as*10, colorS)
legend("Spacecraft Trajectory", "Earth (200 scale)", "Sun Trajectory", "Sun at tf (10 scale)")
axis equal
grid on
hold off


figure
hold on
plotBody([0;0;0], as*10, colorS)
plot3(r_S2E_N(1,:), r_S2E_N(2,:), r_S2E_N(3,:), 'Color', colorE)
plotBody(r_S2E_N(:,end), constants.ae*500, colorE)
plot3(X_2a_true(1,:)+r_S2E_N(1,:), X_2a_true(2,:)+r_S2E_N(2,:), X_2a_true(3,:)+r_S2E_N(3,:))
scatter3(X_2a_true(1,end)+r_S2E_N(1,end), X_2a_true(2,end)+r_S2E_N(2,end), X_2a_true(3,end)+r_S2E_N(3,end), '*')
% legend("Spacecraft Trajectory", "Earth (200 scale)", "Sun Trajectory", "Sun at tf (10 scale)")
axis equal
grid on
hold off


%% part 2 estimate with observations

constants.DSS_latlon(2,2) = -355.749444*pi/180;

iterations = 1;
[X_2a, Xref_2a, dx_2a, P_2a, y_2a, alpha_2a] = CKF_SNC(tspan_2a, data_2a, meas_cov, zeros(3), X0_2a, zeros(n,1), P0_2a, constants);
plotResiduals(tspan_2a, y_2a, alpha_2a, [rngNoise, rngRtNoise], ["2a Pre-Fit Residuals vs. Time", "2a Post-Fit Residuals vs. Time"])
error2a = X_2a-X_2a_true;
plotErrorAndBounds_proj2(tspan_2a, error2a, P_2a, "2a State Error and 3$\sigma$ Bounds vs. Time")

%% part 3
constants.DSS_latlon(2,2) = 355.749444*pi/180;