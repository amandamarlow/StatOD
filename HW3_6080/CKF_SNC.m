function [X, dx, P, y, alpha] = CKF_SNC(t, data, R, Qc, X0, dx0, P0, constants)
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
    % dx(:,1) = dx0;
    P(:,:,1) = P0;
%     STM = eye(n);

    %% Integrate reference trajectory and STM from t(i-1) to t(i)
    [X, STM_t0] = integrateTrajectorySTM_HW3(t, X(:,1), eye(n), constants);
    stationNum = data(1,2);
    [GofXt, H] = GandH_HW3(t(1), X(:,1), stationNum, constants);
    Y = data(1,3:4)';
    y(:,1) = Y - GofXt; % pre-fit residual
    K = P0*H'/(H*P0*H' + R); % kalman gain
    dx(:,1) = dx0 + K*(y(:,1) - H*dx0);
    P(:,:,1) = (eye(n) - K*H)*P0*(eye(n) - K*H)' + K*R*K';
    alpha(:,1) = y(:,1) - H*dx(:,1); % post-fit residual 
    for i = 2:length(t)
        %% Time update
        STM = STM_t0(:,:,i)/STM_t0(:,:,i-1);
        dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        %% SNC
        dt = t(i) - t(i-1);
        Qd = [dt^3/3*Qc, dt^2/2*Qc;
            dt^2/2*Qc, dt*Qc];
        P_ap = STM*P(:,:,i-1)*STM' + Qd;
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == data(:,1));
        if isempty(Y_at_t)
%             dx(:,i) = dx_ap;
%             P(:,:,i) = P_ap;
            y(:,i) = NaN(2,1); % pre-fit residual
            alpha(:,i) = NaN(2,1); % post-fit residual
        else
            for j = Y_at_t
                Y = data(j,3:4)';
                stationNum = data(j,2);
                [GofXt, H] = GandH_HW3(t(i), X(:,i), stationNum, constants);
                y(:,i) = Y - GofXt; % pre-fit residual
                K = P_ap*H' * inv(H*P_ap*H' + R); % kalman gain
                dx_ap = dx_ap + K*(y(:,i) - H*dx_ap);
                P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R*K';
                alpha(:,i) = y(:,i) - H*dx_ap; % post-fit residual 
            end  
        end
        dx(:,i) = dx_ap;
        P(:,:,i) = P_ap;
    end

end

