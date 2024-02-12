function [X, dx, P, y, alpha] = batch(t, Ymat, R, X0, dx0, P0, constants)
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
    STM = eye(n);

    if min(Ymat(:,1) == t(1))
        Y = Ymat(1,3:4)';
        stationNum = Ymat(1,2);
        [GofXt0, H0] = GandH(t(1), X(:,1), stationNum, constants);
        A = P0^-1 + H0'*R(:,:,1)^-1*H0; % Information Matrix
        y(:,1) = Y - GofXt0; % pre-fit residual
        N = H0'*R(:,:,1)^-1*y(:,1);       
    else
        A = P0^-1; % Information Matrix
        N = P0^-1*dx0;
    end
    
    while norm(dx)>1e-12 % is this the right condition?
        for i = 2:length(t)
            %% Integrate reference trajectory and STM from t(i-1) to t(i)
            [X(:,i), STM] = integrateTrajectorySTM(t(i-1), t(i), X(:,i-1), STM, constants);
            %% Time update
%             dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
%             P_ap = STM*P(:,:,i-1)*STM';
            %% Compute observation deviation, observation-state matrix, kalman gain matrix
            Y_at_t = find(t(i) == Ymat(:,1));
            if isempty(Y_at_t)
                y(:,i) = NaN(2,1); % pre-fit residual
                alpha(:,i) = NaN(2,1); % post-fit residual
            else
                for j = Y_at_t
                    Y = Ymat(j,3:4)';
                    stationNum = Ymat(j,2);
                    [GofXt, H] = GandH(t(i), X(:,i), stationNum, constants);
                    if isnan(GofXt(1))
                        y(:,i) = NaN(2,1); % pre-fit residual
                        alpha(:,i) = NaN(2,1); % post-fit residual 
                    else
                        y(:,i) = Y - GofXt; % pre-fit residual
                        A = A + (H*STM)'*R(:,:,j)^-1*H;
                        N = N + (H*STM)'*R(:,:,j)^-1*y;
                    end
                end  
            end
        end
        % solve normal equations
        dx = A^-1*N;
        P = A^-1;
    end
    alpha(:,i) = y(:,1) - H0*STM*dx; % post-fit residual
end

