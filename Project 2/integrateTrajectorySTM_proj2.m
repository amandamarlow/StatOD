function [Xvec, STM_vec] = integrateTrajectorySTM_proj2(tspan, X_t1, STM_t1, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   Detailed explanation goes here

n = length(X_t1);
% tspan = [t1 t2];
S0 = [X_t1; reshape(STM_t1,[],1)];

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
if n > 7
    [~,S] = ode45(@(t,S) proj2ODEwConsts(t,S,constants), tspan, S0, options);
else
    [~,S] = ode45(@(t,S) project2ODE(t,S,constants), tspan, S0, options);
end

% STM
flatSTM = S(:,n+1:end);
STM_vec = reshape(flatSTM',n,n,[]);

% return
Xvec = S(:,1:n)';
% STM_t2 = STM(:,:,end);

end

