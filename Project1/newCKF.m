function [X, dx, P, y, alpha] = newCKF(t, data, R, X0, dx0, P0, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
    t = data(:,1);
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
    for i = 2:length(t)
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [Xvec, STMvec] = integrateTrajectorySTM([t(i-1) t(i)], X(:,i-1), eye(n), constants);
        STM = STMvec(:,:,end);
        X(:,i) = Xvec(:,end);
        %% Time update
        dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        P_ap = STM*P(:,:,i-1)*STM';
        
        Y = data(i,3:4)';
        stationNum = data(i,2);
        [GofXt, H] = GandH(t(i), X(:,i), stationNum, constants);
        y(:,i) = Y - GofXt; % pre-fit residual
        K = P_ap*H'/(H*P_ap*H' + R(:,:,i)); % kalman gain
        dx_ap = dx_ap + K*(y(:,i) - H*dx_ap);
        P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R(:,:,i)*K';
        alpha(:,i) = y(:,i) - H*dx_ap; % post-fit residual 

        dx(:,i) = dx_ap;
        P(:,:,i) = P_ap;
    end

end

