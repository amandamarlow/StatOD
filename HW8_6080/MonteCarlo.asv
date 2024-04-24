function [Xtrue, Xnav] = MonteCarlo(tspan, Xtrue0, Xnav0, Pnav0, targs, S, R, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = 2;
H = [1,0];
B = [0;1];

Xtrue = zeros(n, length(tspan), N);
Xnav = zeros(n, length(tspan), N);
Xtrue(:,1,:) = Xtrue0;
Xnav(:,1,:) = Xnav0;

% w = mvnrnd([0;0],S,N*length(tspan)/2)';
% v = mvnrnd(0,R,N)';
for i = 1:N
    maneuver_num = 1;
    Pnav = Pnav0;
    for k = 2:length(tspan)
        t = tspan(k);
        Phi = STM(tspan(k)-tspan(k-1));
        Xtrue(:,k,i) = Phi*Xtrue(:,k-1,i);
        Xnav(:,k,i) = Phi*Xnav(:,k-1,i);
        Pnav = Phi*Pnav*Phi';
        
        if t==targs(maneuver_num, 1)
            targVec = targs(maneuver_num,2:end);
            dV = waypoint(Xnav(:,k,i), t, targVec);
            % Xtrue(:,k,i) = Xtrue(:,k,i) + B*(dV) + mvnrnd([0;0],S,1)';
            % Xtrue(:,k,i) = Xtrue(:,k,i) + mvnrnd([0,dV],S)';
            Xtrue(:,k,i) = Xtrue(:,k,i) + B*normrnd(dV,S(2,2).^(1/2));
            Xnav(:,k,i) = Xnav(:,k,i) + B*dV;
            Pnav = Pnav + S;
            if maneuver_num < 3
                maneuver_num = maneuver_num+1;
            end

        elseif (mod(t,0.2)==0) && (t~=0)
            % y = H*Xtrue(:,k,i) + mvnrnd(0,R);
            y = H*Xtrue(:,k,i) + normrnd(0,sqrt(R));
            % y = H*Xtrue(:,k,i);
            % residual = y - H*Xnav(:,k,i);
            K = Pnav*H' /(H*Pnav*H' + R); % kalman gain
            % Xnav(:,k,i) = Xnav(:,k,i) + K*(residual - H*Xnav(:,k,i));
            Xnav(:,k,i) = Xnav(:,k,i) + K*(y - H*Xnav(:,k,i));
            Pnav = (eye(n) - K*H)*Pnav*(eye(n) - K*H)' + K*R*K';
        end
        
    end
end

end

