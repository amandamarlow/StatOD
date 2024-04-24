function [C, A] = measUpdate_lincov(X_ap, A_ap, C_ap, B, H, R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = 2;

P_ap = C_ap

K = P_ap*H' / (H*P_ap*H' + R); % kalman gain
A = [eye(n), zeros(n); K*H, eye(n)-K*H];
C = A*C_ap*A' + B*R*B;
end

