function [X, dx, P, y, alpha] = CKF(t, Ymat, R, X0, dx0, P0, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
%     t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(X0);
    X = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    y = zeros(2, length(t)); % pre-fit residual
    alpha = zeros(2, length(t)); % post-fit residual
    %% Initialize
    X(:,1) = X0;
    dx(:,1) = dx0;
    P(:,:,1) = P0;
%     STM = eye(n);
    for i = 2:length(t)
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [X(:,i), STM] = integrateTrajectorySTM(t(i-1), t(i), X(:,i-1), eye(n), constants);
        %% Time update
        dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        P_ap = STM*P(:,:,i-1)*STM';
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == Ymat(:,1));
        if isempty(Y_at_t)
%             dx(:,i) = dx_ap;
%             P(:,:,i) = P_ap;
            y(:,i) = NaN(2,1); % pre-fit residual
            alpha(:,i) = NaN(2,1); % post-fit residual
        else
            for j = Y_at_t
                Y = Ymat(j,3:4)';
                stationNum = Ymat(j,2);
                [GofXt, H] = GandH_HW2(t(i), X(:,i), stationNum, constants);
                if isnan(GofXt(1))
%                     dx(:,i) = dx_ap;
%                     P(:,:,i) = P_ap;
                    y(:,i) = NaN(2,1); % pre-fit residual
                    alpha(:,i) = NaN(2,1); % post-fit residual 
                else
                    y(:,i) = Y - GofXt; % pre-fit residual
                    K = P_ap*H' * inv(H*P_ap*H' + R(:,:,j)); % kalman gain
                    dx_ap = dx_ap + K*(y(:,i) - H*dx_ap);
                    P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R(:,:,j)*K';
                    alpha(:,i) = y(:,i) - H*dx_ap; % post-fit residual 
                end
            end  
        end
        dx(:,i) = dx_ap;
        P(:,:,i) = P_ap;
    end

end

