clear
% clc
close all

addpath('C:\Users\marlo\OneDrive - UCB-O365\Documents\6080 - Statistical Orbit Determination\Project 2\Provided Files')
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
addpath('C:\Users\marlo\MATLAB Drive\6080\StatOD')

% readtable("Project2a_Obs_modified.txt")
data_2 = readmatrix("Project2a_Obs.txt");
data_2(:,end) = [];
data_3 = readmatrix("Project2b_Obs.txt");
data_3(:,end) = [];
% readtable('Project2 Prob2 truth traj 50days.txt')
truthTraj = readmatrix("Project2_Prob2_truth_traj_50days.txt");
truthTraj(:,end) = [];

% setup synamics
t0 = 2456296.25; % JD
constants.t0 = t0;
muS = 132712440017.987; % km^3/s^2 gravitiational parameter of sun
constants.muS = muS;
% constants.muE = 398600.4415; % [km^3/s^2] earth's gravitational parameter
srp_flux = 1357; % W/m^2 at 1 AU
c = 299792458; % m/s
constants.AU = 149597870.7; % km
p_srAU = srp_flux/c*1000; % kN
constants.p_srAU = p_srAU;
a2m = 0.01e-6; % km^2/kg
constants.a2m = a2m;
RSOI = 925000; % km
constants.RSOI = RSOI;
[~, ~, muE] = Ephem(t0,3,'EME2000');

% measurement setup
constants.theta0 = 0;
constants.omegaE = 7.29211585275553e-5; % rad/s
ae = 6378.1363; % km
constants.ae = ae;
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
P0_2 = diag([sigX^2*ones(1,3), sigV^2*ones(1,3), sigCr^2]);
S0_traj = [X0_true; reshape(eye(n),[],1)];

% compare to provided truth trajectory
tspan_traj = truthTraj(:,1);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_traj,S_traj] = ode45(@(t,S) project2ODE(t,S,constants), tspan_traj, S0_traj, options);
r_traj = S_traj(:,1:3)';

% figure
% plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:))
% hold on
% % plotBody([0;0;0], constants.ae*1000)
% plot3(truthTraj(:,2), truthTraj(:,3), truthTraj(:,4), '--')
% axis equal
% grid on
% hold off
% 
% trajError = truthTraj-[t_traj, S_traj];
% final_error = trajError(end, :);
% fprintf("final state errors:")
% disp(final_error(1:8))
% fprintf("final STM errors: \n")
% disp(reshape(final_error(9:end), n,n))

% figure
% subplot(2,1,1)
% semilogy(abs(trajError))
% subplot(2,1,2)
% semilogy(vecnorm(abs(trajError),2,2))

% simulate whole trajectory
as = 696000; % km
tspan_2a = data_2(:,1);
X0_2 = [-274096790.0; -92859240.0; -40199490.0; 32.67; -8.94; -3.88; 1.2]; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
n = length(X0_2);
% S0_2a = [X0_2a; reshape(eye(n),[],1)];
S0_true = [X0_true; reshape(eye(n),[],1)];
[t_2a,S_2a] = ode45(@(t,S) project2ODE(t,S,constants), tspan_2a, S0_true, options);
X_2a_true = S_2a(:,1:n)';

colorS = "#EDB120";
% colorE = "#77AC30";
colorE = [0.4660 0.6740 0.1880];

% r_S2E_N = zeros(3,length(tspan_2a));
% for i = 1:length(tspan_2a)
%     [r_S2E_N(:,i), ~, ~] = Ephem(constants.t0+tspan_2a(i)/86400,3, 'EME2000');
% end
% r_E2S_N = -r_S2E_N;

% figure
% plot3(X_2a_true(1,:), X_2a_true(2,:), X_2a_true(3,:))
% hold on
% plotBody([0;0;0], constants.ae*200, colorE)
% plot3(r_E2S_N(1,:), r_E2S_N(2,:), r_E2S_N(3,:), 'Color', colorS)
% plotBody(-r_S2E_N(:,end), as*10, colorS)
% legend("Spacecraft Trajectory", "Earth (200 scale)", "Sun Trajectory", "Sun at tf (10 scale)")
% axis equal
% grid on
% hold off


% figure
% hold on
% plotBody([0;0;0], as*10, colorS)
% plot3(r_S2E_N(1,:), r_S2E_N(2,:), r_S2E_N(3,:), 'Color', colorE)
% plotBody(r_S2E_N(:,end), constants.ae*500, colorE)
% plot3(X_2a_true(1,:)+r_S2E_N(1,:), X_2a_true(2,:)+r_S2E_N(2,:), X_2a_true(3,:)+r_S2E_N(3,:))
% scatter3(X_2a_true(1,end)+r_S2E_N(1,end), X_2a_true(2,end)+r_S2E_N(2,end), X_2a_true(3,end)+r_S2E_N(3,end), '*')
% % legend("Spacecraft Trajectory", "Earth (200 scale)", "Sun Trajectory", "Sun at tf (10 scale)")
% axis equal
% grid on
% hold off


%% part 2 estimate with observations

constants.DSS_latlon(2,2) = -355.749444*pi/180;
% 
% iterations = 1;
% % [X_2a, Xref_2a, dx_2a, P_2a, y_2a, alpha_2a] = CKF_proj2(tspan_2a, data_2a, meas_cov, zeros(3), X0_2a, zeros(n,1), P0_2a, 1, false, constants);
% % plotResiduals(tspan_2a, y_2a, alpha_2a, [rngNoise, rngRtNoise], ["2a Pre-Fit Residuals vs. Time", "2a Post-Fit Residuals vs. Time"])
% % error2a = X_2a-X_2a_true;
% % plotErrorAndBounds_proj2(tspan_2a, error2a, P_2a, "2a State Error and 3$\sigma$ Bounds vs. Time")
% 
% % [X_2a_smooth, Xref_2a_smooth, dx_2a_smooth, P_2a_smooth, y_2a_smooth, alpha_2a_smooth] = CKF_project2(tspan_2a, data_2a, meas_cov, zeros(3), X0_2a, zeros(n,1), P0_2a, 1, true, constants);
% % plotResiduals(tspan_2a, y_2a_smooth, alpha_2a_smooth, [rngNoise, rngRtNoise], ["2a Pre-Fit Residuals vs. Time", "2a Post-Fit Residuals vs. Time"])
% % error2a = X_2a_smooth-X_2a_true;
% % plotErrorAndBounds_proj2(tspan_2a, error2a_smooth, P_2a_smooth, "2a State Error and 3$\sigma$ Bounds vs. Time")
% 
% Qc = zeros(3);
% % sigma_SNC = 1e-10;
% % Qc = diag(sigma_SNC^2*ones(1,3));
% [X_2a_iterated, Xref_2a_iterated, dx_2a_iterated, P_2a_iterated, y_2a_iterated, alpha_2a_iterated] = CKF_proj2(tspan_2a, data_2, meas_cov, Qc, X0_2, zeros(n,1), P0_2, 15, false, constants);
% plotResiduals(tspan_2a, y_2a_iterated, alpha_2a_iterated, [rngNoise, rngRtNoise], ["2a Pre-Fit Residuals vs. Time", "2a Post-Fit Residuals vs. Time"])
% error2a_iterated = X_2a_iterated-X_2a_true;
% plotErrorAndBounds_proj2(tspan_2a, error2a_iterated, P_2a_iterated, "2a State Error and 3$\sigma$ Bounds vs. Time")
% 
% % notPosDef = zeros(1,length(P_2a_iterated));
% % for i = 1:length(P_2a_iterated)
% %     notPosDef(i) = all(diag(P_2a_iterated(:,:,i))<=0);
% % end
% % sumNotPosDef = sum(notPosDef);
% 
% % sigma_SNC = 0.5e-8;
% % Qc = diag(sigma_SNC^2*ones(1,3));
% % [X_2a_iteratedSmooth, Xref_2a_iteratedSmooth, dx_2a_iteratedSmooth, P_2a_iteratedSmooth, y_2a_iteratedSmooth, alpha_2a_iteratedSmooth] = CKF_proj2(tspan_2a, data_2a, meas_cov, Qc, X0_2a, zeros(n,1), P0_2a, 15, true, constants);
% % plotResiduals(tspan_2a, y_2a_iteratedSmooth, alpha_2a_iteratedSmooth, [rngNoise, rngRtNoise], ["2a Pre-Fit Residuals vs. Time", "2a Post-Fit Residuals vs. Time"])
% % error2a_iteratedSmooth = X_2a_iteratedSmooth-X_2a_true;
% % plotErrorAndBounds_proj2(tspan_2a, error2a_iteratedSmooth, P_2a_iteratedSmooth, "2a State Error and 3$\sigma$ Bounds vs. Time")

%% B-Plane
% t_DCO_2 = tspan_2a(end);
% X_DCO_2 = X_2a_iterated(:,end);
% % X_DCO_2 = X_2a_iteratedSmooth(:,end);

% trueBdotR = 14970.824; % km
% trueBdotT = 9796.737; % km
% % trueBdotVec = [trueBdotR;trueBdotT];
% trueBdotVec = [trueBdotT;trueBdotR];
% 
% [t23RSOI, Xvec_to3RSOI, STM_fromDCO] = integrateTo3RSOI([t_DCO_2, t_DCO_2*100], X_DCO_2, eye(n), constants);
% r3SOI_N = Xvec_to3RSOI(1:3,end);
% v3SOI_N = Xvec_to3RSOI(4:6,end);
% P3SOI_N = STM_fromDCO(:,:,end)*P_2a_iterated(:,:,end)*STM_fromDCO(:,:,end)';
% % P3SOI_N = STM_fromDCO(:,:,end)*P_2a_iteratedSmooth(:,:,end)*STM_fromDCO(:,:,end)';
% P3SOI_N = P3SOI_N(1:6,1:6,end);
% varyShat = true;
% [BdotVec, B_STM] = Bplane(r3SOI_N, v3SOI_N, muE, varyShat);
% % [BdotVec, B_STM] = Bplane(r3SOI_N, v3SOI_N, muE, varyShat);
% P_B = B_STM*P3SOI_N*B_STM';
% 
% fprintf("B-Plane Error: \n")
% disp(BdotVec-trueBdotVec)
% 
% figure
% hold on
% plot1 = plotBplaneCovEllipse(P_B, BdotVec);
% % error_ellipse(P_B, [BdotT; BdotR], 'conf',0.63)
% scatter(BdotVec(1),BdotVec(2), '*')
% scatter(trueBdotVec(1),trueBdotVec(2), '*')
% legend("$3\sigma$ Bounds", "Estimated Target", "True Target", 'Interpreter', 'latex', 'Location','best')
% axis equal
% xlabel("B$\cdot$T", Interpreter="latex")
% ylabel("B$\cdot$R", Interpreter="latex")
% 
% t_DCOtarg_vec = (50:50:200)*24*60^2;
% t_DCO_vec = zeros(1,length(t_DCOtarg_vec));
% X_DCO_vec = zeros(n,length(t_DCO_vec));
% BdotVec_vec = zeros(2,length(t_DCO_vec));
% P_DCO_vec = zeros(n,n,length(t_DCO_vec));
% P_3SOI_vec = zeros(6,6,length(t_DCO_vec));
% P_B_vec = zeros(2,2,length(t_DCO_vec));
% Qc = zeros(3);
% % sigma_SNC = 0.5e-8;
% % Qc = diag(sigma_SNC^2*ones(1,3));
% smooth = false;
% iterations = 15;
% for i = 1:length(t_DCO_vec)
%     tspan_temp = tspan_2a(tspan_2a<=t_DCOtarg_vec(i));
%     t_DCO_vec(i) = tspan_temp(end);
%     [X_looped, ~, ~, P_looped, ~, ~, ~] = CKF_proj2(tspan_temp, data_2, meas_cov, Qc, X0_2, zeros(n,1), P0_2, iterations, smooth, constants);
%     X_DCO_vec(:,i) = X_looped(:,end);
%     P_DCO_vec(:,:,i) = P_looped(:,:,end);
%     [t23RSOI, Xvec_to3RSOI, STM_fromDCO] = integrateTo3RSOI([t_DCO_vec(i), t_DCO_vec(i)*100], X_DCO_vec(:,i), eye(n), constants);
%     r3SOI_N(:,i) = Xvec_to3RSOI(1:3,end);
%     v3SOI_N(:,i) = Xvec_to3RSOI(4:6,end);
%     P_3SOI_temp = STM_fromDCO(:,:,end)*P_DCO_vec(:,:,i)*STM_fromDCO(:,:,end)';
%     P_3SOI_vec(:,:,i) = P_3SOI_temp(1:6,1:6);
% end
% 
% %% B-Plane Plotting
% BdotVec_fixS = zeros(2,length(t_DCO_vec));
% P_B_fixS = zeros(2,2,length(t_DCO_vec));
% BdotVec_varS = zeros(2,length(t_DCO_vec));
% P_B_varS = zeros(2,2,length(t_DCO_vec));
% for i = 1:length(t_DCO_vec)
%     [BdotVec_fixS(:,i), B_STM_fixLoop] = Bplane(r3SOI_N(:,i), v3SOI_N(:,i), muE, false);
%     [BdotVec_varS(:,i), B_STM_varLoop] = Bplane(r3SOI_N(:,i), v3SOI_N(:,i), muE, true);
%     P_B_fixS(:,:,i) = B_STM_fixLoop*P_3SOI_vec(:,:,i)*B_STM_fixLoop';
%     P_B_varS(:,:,i) = B_STM_varLoop*P_3SOI_vec(:,:,i)*B_STM_varLoop';
% end
% legendStrings = ["50 days", "", "100 days", "", "150 days", "", "200 days", ""];
% plotBplane(BdotVec_fixS, P_B_fixS, legendStrings, "B-Plane (Fixed $\hat{S}$) - Question 2h", trueBdotVec, constants.ae)
% plotBplane(BdotVec_varS, P_B_varS, legendStrings, "B-Plane (Variable $\hat{S}$) - Question 2h", trueBdotVec, constants.ae)
% plotBplane(BdotVec_fixS, P_B_fixS, legendStrings, "B-Plane (Fixed $\hat{S}$) - Question 2h", trueBdotVec, 0)
% plotBplane(BdotVec_varS, P_B_varS, legendStrings, "B-Plane (Variable $\hat{S}$) - Question 2h", trueBdotVec, 0)

%% part 3
constants.DSS_latlon(2,2) = 355.749444*pi/180;
tspan_3 = data_3(:,1);
X0_3 = [-274096770.76544; -92859266.4499061; -40199493.6677441; 32.6704564599943; -8.93838913761049; -3.87881914050316; 1.0]; % X km, Y km, Z km, VX km/sec, VY km/sec, VZ km/sec, CR
P0_3 = P0_2;

% rngNoise = 5e-1; % km
% rngRtNoise = 0.5e-4; % km
meas_cov = [rngNoise^2, 0; 0, rngRtNoise^2];

%% regular lkf
% % Qc = zeros(3);
% sigma_SNC = 1e-6; % 1e-7 too little for bad data
% Qc = diag(sigma_SNC^2*ones(1,3));
% [X_3_lkf, Xref_3_lkf, dx_3_lkf, P_3_lkf, y_3_lkf, alpha_3_lkf, cond_lkf] = CKF_proj2(tspan_3, data_3, meas_cov, Qc, X0_3, zeros(n,1), P0_2, 15, false, constants);
% plotResiduals(tspan_3, y_3_lkf, alpha_3_lkf, [rngNoise, rngRtNoise], ["Part 3 Pre-Fit Residuals vs. Time", "Part 3 Post-Fit Residuals vs. Time"])
% X_DCO_lkf = X_3_lkf(:,end);
% P_DCO_lkf = P_3_lkf(:,:,end);
% 
% [~, r3SOI_3, v3SOI_3, P_3SOI] = propogateTo3RSOI(tspan_3(end), X_DCO_lkf, P_DCO_lkf, constants);
% varyShat = false;
% [BdotVec_3, B_STM_3] = Bplane(r3SOI_3, v3SOI_3, muE, varyShat);
% P_B_3 = B_STM_3*P_3SOI*B_STM_3';
% legendStrings = ["3$\sigma$", "Estimated Target"];
% plotBplane(BdotVec_3, P_B_3, legendStrings, "B-Plane (Fixed $\hat{S}$) - Question 2h", [], constants.ae)
% 
% condCheck = log(cond_lkf);
% figure
% plot(condCheck)
% hold on
% yline(15)
% title("Condition Number of $(H*P^-*H' + R)$", 'Interpreter','latex')

%% part 3 with constants
r_stns_E = zeros(3,3);
for i = 1:3
    r_stns_E(:,i) = latlon2ECEF(ae+constants.DSS_alt(i), constants.DSS_latlon(i,1), constants.DSS_latlon(i,2));
end
% X0_3const = [X0_3; muE; muS; p_srAU; a2m; reshape(r_stns_E, [],1)];
% X0_3const = [X0_3; muE; muS; p_srAU; reshape(r_stns_E, [],1)];
% P0_constants = diag([1e-61, 1e-60, 0.00453*10^-2, 1e-10, 1e-16*ones(1,3), 1*ones(1,6)]);
% P0_constants = diag([1e-40, 1e-30, 0.00453*10^-2, 1e-10, 1e-16*ones(1,3), 1*ones(1,6)]);
% P0_constants = diag([1e-61, 1e-60, 0.00453*10^-2, 1e-10*ones(1,3), 1e-2*ones(1,6)]);
% P0_constants = diag([1e-61, 1e-60, (0.00453*10^-2)^2, (1e-10)^2, 1e-16*ones(1,3), 1e-2*ones(1,6)]);
X0_3const = [X0_3; muE; muS; reshape(r_stns_E, [],1)];
P0_constants = diag([1e-40, 1e-30, 1e-16*ones(1,3), 1e-2*ones(1,6)]);
P0_3const = blkdiag(P0_3,P0_constants);
n = 18;
% n = 19;

% Qc = zeros(3);
sigma_SNC = 1e-8; % 1e-6 too large
Qc = diag(sigma_SNC^2*ones(1,3));
[X_3_lkf, Xref_3_lkf, dx_3_lkf, P_3_lkf, y_3_lkf, alpha_3_lkf, cond_lkf] = CKF_proj2(tspan_3, data_3, meas_cov, Qc, X0_3const, zeros(n,1), P0_3const, 15, false, constants);
plotResiduals(tspan_3, y_3_lkf, alpha_3_lkf, [rngNoise, rngRtNoise], ["Part 3 Pre-Fit Residuals vs. Time (Incorporated Constants)", "Part 3 Post-Fit Residuals vs. Time (Incorporated Constants)"])
% X_DCO_lkf = X_3_lkf(:,end);
% P_DCO_lkf = P_3_lkf(:,:,end);

% [~, r3SOI_3, v3SOI_3, P_3SOI] = propogateTo3RSOI(tspan_3(end), X_DCO_lkf, P_DCO_lkf, constants);
% varyShat = false;
% [BdotVec_3, B_STM_3] = Bplane(r3SOI_3, v3SOI_3, muE, varyShat);
% P_B_3 = B_STM_3*P_3SOI*B_STM_3';
% legendStrings = ["3$\sigma$", "Estimated Target"];
% plotBplane(BdotVec_3, P_B_3, legendStrings, "B-Plane (Fixed $\hat{S}$) - Question 2h", [], constants.ae)

condCheck = log(cond_lkf);
figure
plot(condCheck)
hold on
yline(15)
title("Condition Number of $(H*P^-*H' + R)$", 'Interpreter','latex')

%% SRIF

[X_SRIF, P_SRIF, u_SRIF, Xref_SRIF, dx_SRIF, y_SRIF, e_SRIF] = SRIF_procNoise(tspan_3, data_3, X0_3const, zeros(n,1), P0_3const, meas_cov, Qc, constants);
plotResiduals(tspan_3, y_SRIF, e_SRIF, [rngNoise, rngRtNoise], ["SRIF Pre-Fit Residuals vs. Time (Incorporated Constants)", "SRIF e vs. Time (Incorporated Constants)"])

%% 
plotDeltaX(tspan_3, X_3_lkf-X0_3const, "constants")
plotDeltaX(tspan_3, X_SRIF-X0_3const, "constants")