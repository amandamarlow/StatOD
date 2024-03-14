function [dx_smooth, state_smooth, P_smooth] = smoother(t, X_ref, dx, P, P_ap, STM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(X_ref(:,end));
state_smooth = zeros(n, length(t));
dx_smooth = zeros(n, length(t));
P_smooth = zeros(n, n, length(t));

% Initialize
dx_smooth(:,end) = dx(:,end);
state_smooth(:,end) = X_ref(:,end) + dx(:,end);
P_smooth(:,:,end) = P(:,:,end);

for k = length(t)-1:-1:1
    Sk = P(:,:,k)*STM(:,:,k+1).'/P_ap(:,:,k+1);
    dx_smooth(:,k) = dx(:,k) + Sk*(dx_smooth(:,k+1) - STM(:,:,k+1)*dx(:,k));
    state_smooth(:,k) = X_ref(:,k) + dx_smooth(:,k);
    % P_smooth(:,:,k) = P(:,:,k) + Sk*(P(:,:,k+1) - P_ap(:,:,k+1))*Sk.';
    P_smooth(:,:,k) = P(:,:,k) + Sk*(P_smooth(:,:,k+1) - P_ap(:,:,k+1))*Sk.';
end

end

