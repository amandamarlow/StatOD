function [C, Pnav] = LINCOV(t, X0, Dtrue0, Pnav0, H, S, R, targs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = 2;

M = eye(2);
C0 = [Dtrue0, Dtrue0*M'; M*Dtrue0, M*Dtrue0*M'+Pnav0];
C = zeros(n*2,n*2,length(t));
Pnav = zeros(n,n,length(t));
% X = zeros(n,length(t));
% P_ap = Pnav0;
C(:,:,1) = C0;
Pnav(:,:,1) = Pnav0;
Xref = X0;

maneuver_num = 1;

for k = 2:length(t)
    
    % Time Update
    Phi = STM(t(k)-t(k-1));
    Phi_aug = blkdiag(Phi, Phi);
    C(:,:,k) = Phi_aug*C(:,:,k-1)*Phi_aug';
    % P_ap = Phi*Pnav(:,:,k-1);
    Pnav(:,:,k) = Phi*Pnav(:,:,k-1)*Phi';
    Xref = Phi*Xref;
    
    
    % Maneuver
    if t(k)==targs(min([maneuver_num,3]), 1)
        targVec = targs(maneuver_num,2:end);
        [dV, Bhat, G] = waypoint(X0, t(k), targVec);
        M = [eye(n); zeros(n)]; 
        D =  [eye(n), Bhat*G'; zeros(n), eye(n) + Bhat*G'];
        % C(:,:,k) = D*C(:,:,k)*D' + M*S*M';
        C(:,:,k) = D*C(:,:,k)*D' + M*S*M';
        Pnav(:,:,k) = Pnav(:,:,k) + S;
        maneuver_num = maneuver_num+1;

    elseif mod(t(k),0.2) == 0
        % Meas Update
        % K = P_ap*H' / (H*P_ap*H' + R); % kalman gain
        K = Pnav(:,:,k)*H' / (H*Pnav(:,:,k)*H' + R); % kalman gain
        B = [zeros(n,1); K];
        A = [eye(n), zeros(n); K*H, eye(n)-K*H];
        % C(:,:,k) = A*C_ap*A' + B*R*B;
        C(:,:,k) = A*C(:,:,k)*A' + B*R*B';
    end

end

end

