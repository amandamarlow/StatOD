function [P, S, Pxx, Pxc, dx, dxc] = measUpdate_consider(Pcc, dc, P_ap, S_ap, dx_ap, Hx, Hc, y, meas_cov)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(dx_ap);

K = P_ap*Hx' / (Hx*P_ap*Hx' + meas_cov);
P = (eye(n) - K*Hx)*P_ap; % do we need joseph form?
S = (eye(n) - K*Hx)*S_ap - K*Hc;
Pxx = P + S*Pcc*S';
Pxc = S*Pcc;
dx = dx_ap + K*(y - Hx*dx_ap);
dxc = dx + S*dc;
end

