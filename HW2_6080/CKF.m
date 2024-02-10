function [dx, P, y] = CKF(Ymat, R, X0, dx0, P0, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
    t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(X0);
    X = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    y = zeros(2, length(t));
    %% Initialize
    X(:,1) = X0;
    dx(:,1) = dx0;
    P(:,:,1) = P0;
    STM = eye(n);
    for i = 2:length(t)
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [X(:,i), STM] = integrateTrajectorySTM(t(i), t(i-1), X(:,i-1), STM, constants);
        %% Time update
        dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        P_ap = STM*P(:,:,i-1)*STM;
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y = Ymat(i,3:4)';
        stationNum = Ymat(i,2);
        [GofXt, H] = GandH(t(i), X(:,i), stationNum, constants);
        %% Measurement Update
        if isnan(GofXt(1)) % WHAT SHOULD WE ACTUALLY DO IF OUR NL MEASUREMENT PREDICTION SAYS NO STATIONS CAN SEE IT
            dx(:,i) = dx_ap;
            P(:,:,i) = P_ap;
        else
            y(:,i) = Y - GofXt; % pre-fit residual
            K = P_ap*H' / (H*P_ap*H' + R(:,:,i)); % kalman gain
            dx(:,i) = dx_ap + K*(y(i) - H*dx_ap);
            P(:,:,i) = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R(:,:,i)*K';
        end
    end

end

