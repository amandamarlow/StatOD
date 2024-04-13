function [Xvec, STM_t2, theta_t2] = integrateTrajectorySTM_consider(tspan, X_t1, STM_t1, theta_t1, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   Detailed explanation goes here

n = length(X_t1);
q = 1;
% tspan = [t1 t2];
S0 = [X_t1; reshape(STM_t1,[],1); reshape(theta_t1,[],1)];

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
% Xvec = S(:,1:n)';
Xvec = S(end,1:n)';
STM_t2 = STM_vec(:,:,end);
theta_t2 = theta_vec(:,:,end);

end
