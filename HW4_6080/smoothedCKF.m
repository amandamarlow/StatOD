function [dx_smooth, state_smooth, P_smooth] = smoothedCKF(t, data, R, Qc, X0, dx0, P0, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
%     t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(X0);
    X_ref = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    P_ap = zeros(n, n, length(t));
    STM = zeros(n, n, length(t));
    y = zeros(2, length(t)); % pre-fit residual
    alpha = zeros(2, length(t)); % post-fit residual
    %% Initialize
    X_ref(:,1) = X0;
    % dx(:,1) = dx0;
    P(:,:,1) = P0;
    j = 1;
%     STM = eye(n);

    %% Integrate reference trajectory and STM from t(i-1) to t(i)
    [X_ref, STM_t0] = integrateTrajectorySTM_HW3(t, X_ref(:,1), eye(n), constants);
    stationNum = data(1,2);
    [GofXt, H] = GandH_HW3(t(1), X_ref(:,1), stationNum, constants);
    Y = data(1,3:4)';
    y(:,1) = Y - GofXt; % pre-fit residual
    K = P0*H'/(H*P0*H' + R); % kalman gain
    dx(:,1) = dx0 + K*(y(:,1) - H*dx0);
    P(:,:,1) = (eye(n) - K*H)*P0*(eye(n) - K*H)' + K*R*K';
    alpha(:,1) = y(:,1) - H*dx(:,1); % post-fit residual 
    for i = 2:length(t)
        %% Time update
        STM(:,:,i) = STM_t0(:,:,i)/STM_t0(:,:,i-1);
        dx_ap = STM(:,:,i)*dx(:,i-1); % a priori state deviation estimate
        %% SNC
        dt = t(i) - t(i-1);
        if t(i) - data(j,1) > 1000
            Qd = zeros(n);
        else
            Qd = [dt^3/3*Qc, dt^2/2*Qc;
                dt^2/2*Qc, dt*Qc];
        end
        % Qd = [dt^3/3*Qc, dt^2/2*Qc;
        %     dt^2/2*Qc, dt*Qc];
        P_ap(:,:,i) = STM(:,:,i)*P(:,:,i-1)*STM(:,:,i)' + Qd;
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == data(:,1));
        if isempty(Y_at_t)
            dx(:,i) = dx_ap;
            P(:,:,i) = P_ap(:,:,i);
            y(:,i) = NaN(2,1); % pre-fit residual
            alpha(:,i) = NaN(2,1); % post-fit residual
        else
            for j = Y_at_t
                Y = data(j,3:4)';
                stationNum = data(j,2);
                [GofXt, H] = GandH_HW3(t(i), X_ref(:,i), stationNum, constants);
                y(:,i) = Y - GofXt; % pre-fit residual
                K = P_ap(:,:,i)*H' * inv(H*P_ap(:,:,i)*H' + R); % kalman gain
                dx(:,i) = dx_ap + K*(y(:,i) - H*dx_ap);
                P(:,:,i) = (eye(n) - K*H)*P_ap(:,:,i)*(eye(n) - K*H)' + K*R*K';
                alpha(:,i) = y(:,i) - H*dx_ap; % post-fit residual 
            end  
        end
    end

    [dx_smooth, state_smooth, P_smooth] = smoother(t, X_ref, dx, P, P_ap, STM);

end

