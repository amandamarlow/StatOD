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
RMSposition3d = rms(vecnorm(error(1:3,:),2,1), "omitnan")';
RMSvelocity3d = rms(vecnorm(error(4:6,:),2,1), "omitnan")';
RMS = rms(error(1:6,:), 2, "omitnan")';


outsideP = 0;
outsidePxx = 0;
for i = 1:length(error)
    if ~all(error(:,i) < 2*(diag(P(:,:,i)).^(1/2)))
        outsideP = outsideP + 1;
    end
    if ~all(error(:,i) < 2*(diag(Pxx(:,:,i)).^(1/2)))
        outsidePxx = outsidePxx + 1;
    end
end
percentP = outsideP/length(P);
percentPxx = outsidePxx/length(Pxx);

Pc_tf = [Pxx(:,:,end), Pxc(:,:,end); Pxc(:,:,end)', Pcc];
[~, PSI] = integrateTrajectoryPSI_consider(tspan, X0, eye(n), zeros(n,q), constants);
PSI_tf = PSI(:,:,end);
STM_tf = PSI(1:n,1:n,end);
dx0map = STM_tf\dx(:,end);
Pc0map = PSI_tf\Pc_tf/PSI_tf';
% Integrate forward in time
[Xint, PSIint] = integrateTrajectoryPSI_consider(tspan, X0+dx0map, eye(n), zeros(n,q), constants);
Pcint = zeros(n+q,n+q,length(tspan));
for k = 1:length(tspan)
    Pcint(:,:,k) = PSIint(:,:,k)*Pc0map*PSIint(:,:,k)';
end
errorint = Xint - XtrueJ3;
plotErrorAndBounds_2sig(tspan, errorint, Pcint(1:n,1:n,:), "State error and $2 \sigma$ bounds with all measurements")

RMSposition3d_int = rms(vecnorm(errorint(1:3,:),2,1), "omitnan")';
RMSvelocity3d_int = rms(vecnorm(errorint(4:6,:),2,1), "omitnan")';
RMS_int = rms(errorint(1:6,:), 2, "omitnan")';

outsidePint = 0;
for i = 1:length(errorint)
    if ~all(errorint(:,i) < 2*(diag(Pcint(1:n,1:n,i)).^(1/2)))
        outsidePint = outsidePint + 1;
    end
end
percentPint = outsidePint/length(Pcint);