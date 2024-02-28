function [X, P, y, alpha] = EKF_SNC(t, data, R, Qc, X0, P0, warmStart, constants)
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
    tfirst = 2;
    nw = 100; % number of measurements processed as an LKF
    if warmStart
        twarm = t(t<=data(nw,1));
        numtw = length(twarm);
        dx0 = zeros(n,1);
        [Xwarm, dxwarm, P(:,:,1:numtw), y(:,1:numtw), alpha(:,1:numtw)] = CKF_SNC(twarm, data(1:numtw,:), R, Qc, X0, dx0, P0, constants);
        X(:,1:nw) = Xwarm+dxwarm;
        tfirst = numtw+1;
    end
    j = 1;
    for i = tfirst:length(t)
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [Xvec, STMvec] = integrateTrajectorySTM_HW3([t(i-1), t(i)], X(:,i-1), eye(n), constants);
        X(:,i) = Xvec(:,end);
        STM = STMvec(:,:,end);
%         [X(:,i), STM] = integrateTrajectorySTM(t(i-1), t(i), X(:,i-1), eye(n), constants);
        %% Time update
%         dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        %% SNC
        dt = t(i) - t(i-1);
        % if t(i) - data(find(data(:,1) < t(i), 1, 'last'), 1) > 1000
        if t(i) - data(j,1) > 1000
            Qd = zeros(n);
        else
            Qd = [dt^3/3*Qc, dt^2/2*Qc;
                dt^2/2*Qc, dt*Qc];
        end
        P_ap = STM*P(:,:,i-1)*STM' + Qd;
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == data(:,1));
        if isempty(Y_at_t)
            P(:,:,i) = P_ap;
            y(:,i) = NaN(2,1); % pre-fit residual
            alpha(:,i) = NaN(2,1); % post-fit residual
        else
            for j = Y_at_t
                Y = data(j,3:4)';
                stationNum = data(j,2);
                [GofXt, H] = GandH_HW3(t(i), X(:,i), stationNum, constants);
                y(:,i) = Y - GofXt; % pre-fit residual
                K = P_ap*H'/(H*P_ap*H' + R); % kalman gain
                dx = K*y(:,i); % CHECK SIGNS
                P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R*K';
                alpha(:,i) = y(:,i) - H*dx; % post-fit residual 
            end
            P(:,:,i) = P_ap;
            X(:,i) = X(:,i) + dx;
        end
    end

end

