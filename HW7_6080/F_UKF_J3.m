function [X_t2] = F_UKF(t1, t2, X_t1, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   assumes state X_t1 is postion and velocity augmented with process noise
%   and measurement noise that are zero mean

% n = length(X_t1);
n = 6;
tspan = [t1 t2];

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,S] = ode45(@(t,S) UKF_J3_ODE(t,S,n,constants), tspan, X_t1(1:n), options);

% return
X_t2 = zeros(length(X_t1),1);
X_t2(1:n) = S(end,1:n)';

end

