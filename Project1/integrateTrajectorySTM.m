function [X_t2, STM_t2] = integrateTrajectorySTM(t1, t2, X_t1, STM_t1, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   Detailed explanation goes here

n = length(X_t1);
tspan = [t1 t2];
S0 = [X_t1; reshape(STM_t1,[],1)];

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,S] = ode45(@(t,S) dragJ2ODE(t,S,constants), tspan, S0, options);

% STM
flatSTM = S(end,n+1:end);
STM_t2 = reshape(flatSTM',n,n,[]);

% return
X_t2 = S(end,1:n)';
% STM_t2 = STM(:,:,end);

end

