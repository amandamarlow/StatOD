function [X, P, y, alpha] = EKF(t, Ymat, R, X0, P0, constants)
%Extended Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
    n = length(X0);
    X = zeros(n, length(t));
    P = zeros(n, n, length(t));
    y = zeros(2, length(t)); % pre-fit residual
    alpha = zeros(2, length(t)); % post-fit residual
    %% Initialize
    X(:,1) = X0;
    P(:,:,1) = P0;
%     STM = eye(n);
    for i = 2:length(t)
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [X(:,i), STM] = integrateTrajectorySTM(t(i-1), t(i), X(:,i-1), eye(n), constants);
        %% Time update
%         dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        P_ap = STM*P(:,:,i-1)*STM';
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == Ymat(:,1));
        dx = zeros(n,1);
        if isempty(Y_at_t)
            P(:,:,i) = P_ap;
            y(:,i) = NaN(2,1); % pre-fit residual
            alpha(:,i) = NaN(2,1); % post-fit residual
        else
            for j = Y_at_t
                Y = Ymat(j,3:4)';
                stationNum = Ymat(j,2);
                [GofXt, H] = GandH(t(i), X(:,i), stationNum, constants);
                if isnan(GofXt(1))
                    P(:,:,i) = P_ap;
                    y(:,i) = NaN(2,1); % pre-fit residual
                    alpha(:,i) = NaN(2,1); % post-fit residual 
                else
                    y(:,i) = Y - GofXt; % pre-fit residual
                    K = P_ap*H' * inv(H*P_ap*H' + R(:,:,j)); % kalman gain
                    dx = dx + K*(y(:,i)- H*dx); % CHECK SIGNS
                    P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R(:,:,j)*K';
                    alpha(:,i) = y(:,i) - H*dx; % post-fit residual 
                end
            end
            P(:,:,i) = P_ap;
            X(:,i) = X(:,i) + dx;
        end
    end

end

