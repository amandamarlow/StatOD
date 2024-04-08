function [Xvec, PSI] = integrateTrajectoryPSI_consider(tspan, X_t0, STM_t0, theta_t0, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   Detailed explanation goes here

n = length(X_t0);
q = 1;
% tspan = [t1 t2];
S0 = [X_t0; reshape(STM_t0,[],1); reshape(theta_t0,[],1)];

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,S] = ode45(@(t,S) considerODE(t,S,n,q,constants), tspan, S0, options);

% STM
flatSTM = S(:,n+1:end-n*q);
STM_vec = reshape(flatSTM',n,n,[]);

% theta
flatTheta = S(:,end-n*q+1:end);
theta_vec = reshape(flatTheta',n,q,[]);

% return
Xvec = S(:,1:n)';
PSI = [STM_vec, theta_vec; zeros(q,n,length(tspan)), ones(1,1,length(tspan))];
end

