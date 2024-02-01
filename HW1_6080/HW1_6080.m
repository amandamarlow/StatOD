clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
load("HW1_truth_traj_mu_J2_with_STM.txt")
load("HW1_truth_J2_J3_acc.txt")

% given orbit elements w/ respect to ECI
a = 10000; %km
e = 0.001; 
i = 40*pi/180; %rad
OMEGA = 80*pi/180; %rad
omega = 40*pi/180; %rad
nu0 = 0;

% constants
mu = 398600.4415;
ae = 6378.0; % [km] mean equitorial radius of earth
% constants.mu = mu; % [km^3/s^2] earth's gravitational parameter
constants.ae = ae; % [km] mean equitorial radius of earth
% constants.J2 = 1.08263*10^-3; % J2 perturbation
J2 = 1.08269*10^-3; % J2 perturbation
% J3 = -0.0000025323;
J3 = 0;
C = [mu; J2; J3];

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
% tspan = [0 15*T];
tspan = HW1_truth_traj_mu_J2_with_STM(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t1,S1] = ode45(@(t,S) orbitODE(t,S,constants), tspan, S0, options);
plotEarthOrbit(S1(:,1:3)', ae, "Simulated Reference Trajectory")

% simulate perturbed state
delta_S0 = [1;0;0;0;0.01;0];
perturbedS = [S0(1:6)+delta_S0; S0(7:end)];
[t2,S2] = ode45(@(t,S) orbitODE(t,S,constants), tspan, perturbedS, options);
true_deviation = S2(:,1:6)' - S1(:,1:6)';

% propogate perturbation using STM
flatSTM = S1(:,10:end);
STM = reshape(flatSTM',9,9,[]);
STM_deviation = zeros(9,length(t1));
STM_deviation(1:6,1) = delta_S0;
for i = 1:length(t1)-1 
   STM_deviation(:,i+1) = STM(:,:,i)*STM_deviation(:,1);
end
comparing_flat_STM = zeros(length(t1),36);
for j = 1:length(t1)
  comparing_flat_STM(j,:) = reshape(STM(1:6,1:6,j),1,[]);
end


error = S1(:,1:6) - HW1_truth_traj_mu_J2_with_STM(:,2:7);
error = [S1(:,1:6), comparing_flat_STM] - HW1_truth_traj_mu_J2_with_STM(:,2:end);

figure
subplot(6, 1, 1);
sgtitle('Difference in deviation vectors vs time');
plot(t2, STM_deviation(1, :)-true_deviation(1,:));
% label
xlabel('time [s]');
ylabel("$\delta x_{STM} - \delta x_{true}$", 'Interpreter', 'latex');
% subplot
subplot(6, 1, 2);
plot(t2, STM_deviation(2, :)-true_deviation(2,:));
% label
xlabel('time [s]');
ylabel("$\delta y_{STM} - \delta y_{true}$", 'Interpreter', 'latex');
% subplot
subplot(6, 1, 3);
plot(t2, STM_deviation(3, :)-true_deviation(3,:));
% label
xlabel('time [s]');
ylabel("$\delta z_{STM} - \delta z_{true}$", 'Interpreter', 'latex');
% subplot
subplot(6, 1, 4);
plot(t2, STM_deviation(4, :)-true_deviation(4,:));
% label
xlabel('time [s]');
ylabel("$\delta \dot{x}_{STM} - \delta \dot{x}_{true}$", 'Interpreter', 'latex');
% subplot
subplot(6, 1, 5);
plot(t2, STM_deviation(5, :)-true_deviation(5,:));
% label
xlabel('time [s]');
ylabel("$\delta \dot{y}_{STM} - \delta \dot{y}_{true}$", 'Interpreter', 'latex');
% subplot
subplot(6, 1, 6);
plot(t2, STM_deviation(6, :)-true_deviation(6,:));
% label
xlabel('time [s]');
ylabel("$\delta \dot{z}_{STM} - \delta \dot{z}_{true}$", 'Interpreter', 'latex');