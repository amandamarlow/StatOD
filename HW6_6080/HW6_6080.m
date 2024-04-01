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

n = 6;
q = 1;
Pxx0 = eye(n)*1e4;
c = 0;
J3 = -0.0000025323;
deltac = J3;
Pcc = J3^2;

%% implement consider covariance filter
dx0 = zeros(n,1);

[X, Xref, dx, P, Pxx, Pxc, y, alpha] = seqConsiderCov(tspan, YJ3_simulated, meas_cov, X0, dx0, P0, Pxx0, deltac, Pcc, constants); % error in Pxx
error = X - XtrueJ3;
plotErrorAndBounds_HW6(tspan, error, P, Pxx, "state error and $2 \sigma$ bounds")

Pc_tf = [Pxx(:,:,end), Pxc(:,:,end); Pxc(:,:,end)', Pcc];
[~, PSI] = integrateTrajectoryPSI_consider(tspan, X0(:,1), eye(n), zeros(n,q), constants);
PSI_tf = PSI(:,:,end);
STM_tf = PSI(1:n,1:n,end);
dx0map = STM_tf\dx(:,end);
Pc0map = PSI_tf\Pc_tf/PSI_tf';
