function [dV, Bhat, G] = waypoint(X0, t, targVec)
%WAYPOINT Summary of this function goes here
%   targVec = [position target, time target]

r_targ = targVec(1);
t_targ = targVec(2);
dt_targ = t_targ - t;

r0 = X0(1);
v0 = X0(2);
Phi = STM(dt_targ);
v_des = Phi(1,2) \ (r_targ - Phi(1,1)*r0);
dV = v_des - v0;

G = -[Phi(1,2)\Phi(1,1); 1];
Bhat = [0;1];
% BG = Bhat*G;
end


