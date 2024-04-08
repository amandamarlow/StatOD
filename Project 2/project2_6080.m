clear
clc
close all

addpath('C:\Users\marlo\OneDrive - UCB-O365\Documents\6080 - Statistical Orbit Determination\Project 2\Provided Files')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')

% readtable("Project2a_Obs_modified.txt")
data_2a = readmatrix("Project2a_Obs.txt");
% readtable('Project2 Prob2 truth traj 50days.txt')
truthTraj = readmatrix("Project2_Prob2_truth_traj_50days.txt");
truthTraj(:,end) = [];

% setup synamics
constants.t0 = 2456296.25; % JD
constants.muS = 132712440017.987; % km^3/s^2 gravitiational parameter of sun
constants.muE = 398600.4415; % [km^3/s^2] earth's gravitational parameter
srp_flux = 1357; % W/m^2 at 1 AU
c = 299792458; % m/s
constants.AU = 149597870.7; % km
constants.p_srAU = srp_flux/c*1000; % kN
constants.a2m = 0.01e-6; % km^2/kg


% measurement setup
theta0 = 0;
omegaE_N = [0; 0; 7.29211585275553e-5]; % rad/s
constants.ae = 6378.1363; % km
DSS_34 = [35.398333*pi/180; 148.981944*pi/180; 0.691750]; % [lat rad; lon rad; altitude km]
DSS_65 = [40.427222*pi/180; 355.749444*pi/180; 0.834539]; % [lat rad; lon rad; altitude km] NEGATIVE ONLY IN PART 2
DSS_13 = [35.247164*pi/180; 243.205*pi/180; 1.07114904]; % [lat rad; lon rad; altitude km]

rngNoise = 5e-3; % km
rngeRtNoise = 0.5e-6; % km
meas_cov = [rngeRtNoise, 0; 0, rngeRtNoise];

%% Estimate the State with a Known Target and Models

% X0 = [-274096790.0; -92859240.0; -40199490.0; 32.67; -8.94; -3.88; 1.0]; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
X0 = truthTraj(1,2:8)'; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
n = length(X0);
sigX = 100; % km
sigV = 0.1; % km/s
sigCr = 0.1;
P0 = diag([sigX^2*ones(1,3), sigV^2*ones(1,3), sigCr^2]);
S0 = [X0; reshape(eye(n),[],1)];

tspan = truthTraj(:,1);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,S] = ode45(@(t,S) project2ODE(t,S,constants), tspan, S0, options);
r_traj = S(:,1:3)';

figure
plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:))
hold on
% plotBody([0;0;0], constants.ae*1000)
plot3(truthTraj(:,2), truthTraj(:,3), truthTraj(:,4), '--')
axis equal
grid on
hold off

trajError = truthTraj-[t, S];