clear
clc
close all

addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD\HW2_6080')
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

n = 6; % length of state vector

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
load("Y_simulated.mat")
tspanIdeal = Y_simulated(:,1);
load("YJ3_simulated.mat")
% tspan = YJ3_simulated(:,1);

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tideal,Sideal] = ode45(@(t,S) orbitODE(t,S,constants), tspanIdeal, S0, options);
[t,S] = ode45(@(t,S) orbitJ3ODE(t,S,constants), tspan, S0, options);
Xtrue = S(:,1:n)';
t = t';


%% Check Sigma Range

X0 = S0(1:n);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];

warmStart = 1;

% exponent = -18:-5; % km/s^2
exponent = -10:-2; % km/s^2
% exponent = -7:-5; % km/s^2

sigsRMSposition_CKF_SNC = zeros(1,length(exponent));
sigsRMSvelocity_CKF_SNC = zeros(1,length(exponent));
sigsRMSresiduals_CKF_SNC = zeros(2,length(exponent));
sigsRMSposition_EKF_SNC = zeros(1,length(exponent));
sigsRMSvelocity_EKF_SNC = zeros(1,length(exponent));
sigsRMSresiduals_EKF_SNC = zeros(2,length(exponent));
sigmas_SNC = zeros(1,length(exponent));
for i = 1:length(exponent)
    sigmas_SNC(i) = 10^exponent(i);
    Qc = diag(sigmas_SNC(i)^2*ones(1,3));

    [tempX_CKF_SNC, tempdx_CKF_SNC, ~, tempy_CKF_SNC, tempalpha_CKF_SNC] = CKF_SNC(t, YJ3_simulated, R, Qc, X0, zeros(6,1), P0, constants);
    temperror_CKF_SNC = (tempX_CKF_SNC + tempdx_CKF_SNC) - Xtrue;
    sigsRMSposition_CKF_SNC(:,i) = rms(vecnorm(temperror_CKF_SNC(1:3,:),2,1), "omitnan");
    sigsRMSvelocity_CKF_SNC(:,i) = rms(vecnorm(temperror_CKF_SNC(4:6,:),2,1), "omitnan");
    sigsRMSresiduals_CKF_SNC(:,i) = rms(tempalpha_CKF_SNC, 2, "omitnan");

    % [tempX_EKF_SNC, ~, tempy_EKF_SNC, tempalpha_EKF_SNC] = EKF_SNC(t, YJ3_simulated, R, Qc, X0, P0, warmStart, constants);
    % temperror_EKF_SNC = tempX_EKF_SNC - Xtrue;
    % sigsRMSposition_EKF_SNC(:,i) = rms(vecnorm(temperror_EKF_SNC(1:3,:),2,1), "omitnan");
    % sigsRMSvelocity_EKF_SNC(:,i) = rms(vecnorm(temperror_EKF_SNC(4:6,:),2,1), "omitnan");
    % sigsRMSresiduals_EKF_SNC(:,i) = rms(tempalpha_EKF_SNC, 2, "omitnan");
end

% a_DMC = zeros(3,1);
% X0 = S0(1:n);
% X0_DMC = [X0; a_DMC];
% P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
% noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
% R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
% tau = T/30;
% warmStart = 1;
% exponent = -16:-8; % km/s^2
% sigsRMSposition_CKF_DMC = zeros(1,length(exponent));
% sigsRMSvelocity_CKF_DMC = zeros(1,length(exponent));
% sigsRMSresiduals_CKF_DMC = zeros(2,length(exponent));
% sigsRMSposition_EKF_DMC = zeros(1,length(exponent));
% sigsRMSvelocity_EKF_DMC = zeros(1,length(exponent));
% sigsRMSresiduals_EKF_DMC = zeros(2,length(exponent));
% sigmas_DMC = zeros(1,length(exponent));
% for i = 1:length(exponent)
%     sigmas_DMC(i) = 10^exponent(i);
%     Qc = diag(sigmas_DMC(i)^2*ones(1,3));
% 
%     P0_aDMC = diag(ones(3,1) * 3*(sigmas_DMC(i))^2*tau/2);
%     P0_DMC = blkdiag(P0, P0_aDMC);
% 
%     [tempX_CKF_DMC, tempdx_CKF_DMC, ~, tempy_CKF_DMC, tempalpha_CKF_DMC] = CKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, zeros(9,1), P0_DMC, tau, constants);
%     temperror_CKF_DMC = (tempX_CKF_DMC(1:6,:)+ tempdx_CKF_DMC(1:6,:)) - Xtrue;
%     sigsRMSposition_CKF_DMC(:,i) = rms(vecnorm(temperror_CKF_DMC(1:3,:),2,1), "omitnan");
%     sigsRMSvelocity_CKF_DMC(:,i) = rms(vecnorm(temperror_CKF_DMC(4:6,:),2,1), "omitnan");
%     sigsRMSresiduals_CKF_DMC(:,i) = rms(tempalpha_CKF_DMC, 2, "omitnan");
% 
%     [tempX_EKF_DMC, ~, tempy_EKF_DMC, tempalpha_EKF_DMC] = EKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, P0_DMC, tau, warmStart, constants);
%     temperror_EKF_DMC = tempX_EKF_DMC(1:6,:) - Xtrue;
%     sigsRMSposition_EKF_DMC(:,i) = rms(vecnorm(temperror_EKF_DMC(1:3,:),2,1), "omitnan");
%     sigsRMSvelocity_EKF_DMC(:,i) = rms(vecnorm(temperror_EKF_DMC(4:6,:),2,1), "omitnan");
%     sigsRMSresiduals_EKF_DMC(:,i) = rms(tempalpha_EKF_DMC, 2, "omitnan");
% end

%% Plot Sigma Range

plotRMSsigmas(sigmas_SNC, sigsRMSposition_CKF_SNC, sigsRMSvelocity_CKF_SNC, sigsRMSresiduals_CKF_SNC, "CKF-SNC RMS Values")
% plotRMSsigmas(sigmas_SNC, sigsRMSposition_EKF_SNC, sigsRMSvelocity_EKF_SNC, sigsRMSresiduals_EKF_SNC, "EKF-SNC RMS Values")
% plotRMSsigmas(sigmas_DMC, sigsRMSposition_CKF_DMC, sigsRMSvelocity_CKF_DMC, sigsRMSresiduals_CKF_DMC, "CKF-DMC RMS Values")
% plotRMSsigmas(sigmas_DMC, sigsRMSposition_EKF_DMC, sigsRMSvelocity_EKF_DMC, sigsRMSresiduals_EKF_DMC, "EKF-DMC RMS Values")



%% Implement Filters with ideal sigma
n = 6;
X0 = S0(1:n);
P0 = blkdiag(eye(3), (1e-3)^2*eye(3));
noise_sd = [1e-3; 1e-6]; % std deviation of the range and range rate noise
R = [noise_sd(1)^2, 0; 0, noise_sd(2)^2];
warmStart = 1;


% sigma_SNC = 10e-6;
sigma_SNC = 1e-6;
Qc = diag(sigma_SNC^2*ones(1,3));

[X_CKF_SNC, dx_CKF_SNC, P_CKF_SNC, y_CKF_SNC, alpha_CKF_SNC] = CKF_SNC(t, YJ3_simulated, R, Qc, X0, zeros(n,1), P0, constants);
error_CKF_SNC = (X_CKF_SNC + dx_CKF_SNC) - Xtrue;

% % [X_EKF_SNC, P_EKF_SNC, y_EKF_SNC, alpha_EKF_SNC] = EKF_SNC(t, YJ3_simulated, R, Qc, X0, P0, warmStart, constants);
% % error_EKF_SNC = X_EKF_SNC - Xtrue;
% 
% n = 9;
% tau = T/30;
% % tau = T/60;
% sigma_DMC = 1e-9;
% Qc = diag(sigma_DMC^2*ones(1,3));
% P0_aDMC = diag(ones(3,1) * 3*(sigma_DMC)^2*tau/2);
% a_DMC = zeros(3,1);
% X0_DMC = [X0; a_DMC];
% P0_DMC = blkdiag(P0, P0_aDMC);
% 
% [X_CKF_DMC, dx_CKF_DMC, P_CKF_DMC, y_CKF_DMC, alpha_CKF_DMC] = CKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, zeros(n,1), P0_DMC, tau, constants);
% error_CKF_DMC = (X_CKF_DMC(1:6,:)+ dx_CKF_DMC(1:6,:)) - Xtrue;
% 
% [X_EKF_DMC, P_EKF_DMC, y_EKF_DMC, alpha_EKF_DMC] = EKF_DMC(t, YJ3_simulated, R, Qc, X0_DMC, P0_DMC, tau, warmStart, constants);
% error_EKF_DMC = X_EKF_DMC(1:6,:) - Xtrue;

%% RMS

% RMSposition_CKF_SNC = rms(vecnorm(error_CKF_SNC(1:3,:),2,1), "omitnan");
% RMSvelocity_CKF_SNC = rms(vecnorm(error_CKF_SNC(4:6,:),2,1), "omitnan");
% RMSresiduals_CKF_SNC = rms(alpha_CKF_SNC, 2, "omitnan");
% 
% 
% RMSposition_EKF_SNC = rms(vecnorm(error_EKF_SNC(1:3,:),2,1), "omitnan");
% RMSvelocity_EKF_SNC = rms(vecnorm(error_EKF_SNC(4:6,:),2,1), "omitnan");
% RMSresiduals_EKF_SNC = rms(alpha_EKF_SNC, 2, "omitnan");


% RMSposition_CKF_DMC = rms(vecnorm(error_CKF_DMC(1:3,:),2,1), "omitnan");
% RMSvelocity_CKF_DMC = rms(vecnorm(error_CKF_DMC(4:6,:),2,1), "omitnan");
% RMSresiduals_CKF_DMC = rms(alpha_CKF_DMC, 2, "omitnan");
% 
% 
% RMSposition_EKF_DMC = rms(vecnorm(error_EKF_DMC(1:3,:),2,1), "omitnan");
% RMSvelocity_EKF_DMC = rms(vecnorm(error_EKF_DMC(4:6,:),2,1), "omitnan");
% RMSresiduals_EKF_DMC = rms(alpha_EKF_DMC, 2, "omitnan");

%% Plot - SNC
titles = ["CKF-SNC Pre-Fit Residuals vs. Time", "CKF-SNC Post-Fit Residuals vs. Time"];
plotResiduals_HW3(t, y_CKF_SNC, alpha_CKF_SNC, noise_sd, titles)
plotErrorAndBounds_HW3(t, error_CKF_SNC, P_CKF_SNC, "CKF-SNC Error vs Time")

% titles = ["EKF Pre-Fit Residuals vs. Time (J3)", "EKF-SNC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_EKF_SNC, alpha_EKF_SNC, noise_sd, titles)
% plotErrorAndBounds_HW3(t, error_EKF_SNC, P_EKF_SNC, "EKF-SNC Error vs Time")

%% Plot - DMC

% titles = ["CKF-DMC Pre-Fit Residuals vs. Time", "CKF-DMC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_CKF_DMC, alpha_CKF_DMC, noise_sd, titles)
% plotErrorAndBounds_HW3(t, error_CKF_DMC, P_CKF_DMC, "CKF-DMC Error vs Time")
% 
% titles = ["EKF-DMC Pre-Fit Residuals vs. Time", "EKF-DMC Post-Fit Residuals vs. Time"];
% plotResiduals_HW3(t, y_EKF_DMC, alpha_EKF_DMC, noise_sd, titles)
% plotErrorAndBounds_HW3(t, error_EKF_DMC, P_EKF_DMC, "EKF-DMC Error vs Time")
% 
% %% a error
% sig_aCKF = zeros(3, length(t));
% aJ3 = zeros(3, length(t));
% P_aCKF = zeros(3,3,length(t));
% for i = 1:length(t)
%     S = X_CKF_DMC(:,i)+dx_CKF_DMC(:,i);
%     x = S(1);
%     y = S(2);
%     z = S(3);
%     r_N = S(1:3);
%     r = norm(r_N);
%     v_N = S(4:6);
% 
%     P_aCKF(:,:,i) = P_CKF_DMC(7:9,7:9,i);
%     aJ3(:,i) = 1/2*mu/r^2*(ae/r)^3*J3 * [5*(7*(z/r)^3 - 3*(z/r))*x/r; 5*(7*(z/r)^3 - 3*(z/r))*y/r; 3*(1 - 10*(z/r)^2 + 35/5*(z/r)^4)];
% end
% error_aCKF = X_CKF_DMC(7:9,:)+dx_CKF_DMC(7:9,:) - aJ3;
% plotJ3aError(t, error_aCKF, P_aCKF, "CKF-DMC Acceleration Error vs Time")
% 
% sig_aEKF = zeros(3, length(t));
% aJ3 = zeros(3, length(t));
% P_aEKF = zeros(3,3,length(t));
% for i = 1:length(t)
%     S = X_EKF_DMC(:,i);
%     x = S(1);
%     y = S(2);
%     z = S(3);
%     r_N = S(1:3);
%     r = norm(r_N);
%     v_N = S(4:6);
% 
%     P_aEKF(:,:,i) = P_EKF_DMC(7:9,7:9,i);
%     aJ3(:,i) = 1/2*mu/r^2*(ae/r)^3*J3 * [5*(7*(z/r)^3 - 3*(z/r))*x/r; 5*(7*(z/r)^3 - 3*(z/r))*y/r; 3*(1 - 10*(z/r)^2 + 35/5*(z/r)^4)];
% end
% error_aEKF = X_EKF_DMC(7:9,:) - aJ3;
% plotJ3aError(t, error_aEKF, P_aCKF, "EKF-DMC Acceleration Error vs Time")
% 
% %% RIC Frame
% 
% % sigma_SNC = 10e-6;
% % Qc_RTN = diag(sigma_SNC^2*ones(1,3));
% % RN = ECI2RTN(r0_ECI, v0_ECI);
% % NR = RN';
% % Qctilde = NR*Qc_RTN*RN;
% % 
% % [X_CKF_SNC, dx_CKF_SNC, P_CKF_SNC, y_CKF_SNC, alpha_CKF_SNC] = CKF_SNC(t, YJ3_simulated, R, Qctilde, X0, zeros(n,1), P0, constants);
% % error_CKF_SNC = (X_CKF_SNC + dx_CKF_SNC) - Xtrue;
% % 
% % [X_EKF_SNC, P_EKF_SNC, y_EKF_SNC, alpha_EKF_SNC] = EKF_SNC(t, YJ3_simulated, R, Qctilde, X0, P0, warmStart, constants);
% % error_EKF_SNC = X_EKF_SNC - Xtrue;
% % 
% % RMSposition_CKF_SNC_RTN = rms(vecnorm(error_CKF_SNC(1:3,:),2,1), "omitnan");
% % RMSvelocity_CKF_SNC_RTN = rms(vecnorm(error_CKF_SNC(4:6,:),2,1), "omitnan");
% % RMSresiduals_CKF_SNC_RTN = rms(alpha_CKF_SNC, 2, "omitnan");
% % 
% % RMSposition_EKF_SNC_RTN = rms(vecnorm(error_EKF_SNC(1:3,:),2,1), "omitnan");
% % RMSvelocity_EKF_SNC_RTN = rms(vecnorm(error_EKF_SNC(4:6,:),2,1), "omitnan");
% % RMSresiduals_EKF_SNC_RTN = rms(alpha_EKF_SNC, 2, "omitnan");
