clear
clc
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

n = 6;
% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,S] = ode45(@(t,S) orbitODE(t,S,constants), tspan, S0, options);
[tJ3,SJ3] = ode45(@(t,S) orbitJ3ODE(t,S,constants), tspan, S0, options);
Xtrue = S(:,1:n)';
XtrueJ3 = SJ3(:,1:n)';
t = t';

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
range_noise = mvnrnd(0,(noise_sd(1))^2, length(range_observations));
rangeRate_noise = mvnrnd(0,(noise_sd(2))^2, length(rangeRate_observations));
Y_simulated = [range_observations(:,1:2), range_observations(:,3) + range_noise, rangeRate_observations(:,3) + rangeRate_noise]; % [t(s), station#, range measurement, range rate measurement]

% simulate measurements with J3
[J3range_observations, J3rangeRate_observations, J3elevations_byStation, J3elevations_all] = simMeas(tJ3, SJ3', [r_s1_E,r_s2_E,r_s3_E], constants);
YJ3_simulated = [J3range_observations(:,1:2), J3range_observations(:,3) + range_noise(1:length(J3range_observations)), J3rangeRate_observations(:,3) + rangeRate_noise(1:length(J3range_observations))]; % [t(s), station#, range measurement, range rate measurement]

%% save data
save('constants', 'constants')
save('Y_simulated', 'Y_simulated')
save('YJ3_simulated', 'YJ3_simulated')
save('Xtrue', 'Xtrue')
save('XtrueJ3', 'XtrueJ3')

% Initial Condition
X0 = S0(1:n);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
meas_cov = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
save('initialConds', 'X0', 'P0', 'noise_sd', 'meas_cov', 'tspan')
