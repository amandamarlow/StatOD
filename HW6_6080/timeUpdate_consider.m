function [P_ap, S_ap, dx_ap, dxc_ap, Pxx_ap, Pxc_ap] = timeUpdate_consider(Pcc, dc, P0, dx0, S0, STM, theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

P_ap = STM*P0*STM'; % covariance of state without consider
S_ap = STM*S0 + theta;
dx_ap = STM*dx0;
dxc_ap = dx_ap + S_ap*dc; 
Pxx_ap = P_ap + S_ap*Pcc*S_ap'; % Covariance of state with consider portion
Pxc_ap = S_ap*Pcc;
end

