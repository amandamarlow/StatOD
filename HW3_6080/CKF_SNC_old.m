function [X, dx, P, y, alpha] = CKF_SNC_old(data, R, Qc, X0, dx0, P0, iterations, constants)
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
%     dx_ap = dx0;
    P(:,:,1) = P0;
%     RMSresidual_post = zeros(2,1);
    RMSresidual_pre = zeros(2,1);
    
    for q = 1:iterations
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [X, STM_t0] = integrateTrajectorySTM_HW3(t, X(:,1), eye(n), constants);
        stationNum = data(1,2);
        [GofXt, H] = GandH_HW3(t(1), X(:,1), stationNum, constants);
        Y = data(1,3:4)';
        y(:,1) = Y - GofXt; % pre-fit residual
%         K = P0*H'/(H*P0*H' + R); % kalman gain
%         dx(:,1) = dx_ap + K*(y(:,1) - H*dx_ap);
%         P(:,:,1) = (eye(n) - K*H)*P0*(eye(n) - K*H)' + K*R*K';
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
            %% Measurement Update
            Y = data(i,3:4)';
            stationNum = data(i,2);
            [GofXt, H] = GandH_HW3(t(i), X(:,i), stationNum, constants);
            y(:,i) = Y - GofXt; % pre-fit residual
            K = P_ap*H'/(H*P_ap*H' + R); % kalman gain
            dx_ap = dx_ap + K*(y(:,i) - H*dx_ap);
            P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R*K';
            alpha(:,i) = y(:,i) - H*dx_ap; % post-fit residual 

            dx(:,i) = dx_ap;
            P(:,:,i) = P_ap;
        end
%         RMSresidual_qminus = RMSresidual_post;
        RMSresidual_qminus = RMSresidual_pre;
        RMSresidual_pre = sqrt(1/length(t)*sum(y.^2, 2));
%         RMSresidual_post = sqrt(1/length(t)*sum(alpha.^2, 2));
%         if all(abs(RMSresidual_post-RMSresidual_qminus) < 0.01*RMSresidual_post) && q>1
        if all(abs(RMSresidual_pre-RMSresidual_qminus) < 0.01*RMSresidual_pre) && q>1
            fprintf("CKF Converged in %d iterations \n", q)
            break
        elseif q == iterations && iterations>1
            fprintf("CKF Failed to converged in %d iterations \n", q)
            break
        elseif q == iterations
            break
        end
        STM = STM_t0(:,:,end);
        diffX0 = STM\dx(:,end);
%         dx(:,1) = dx(:,1) - diffX0;
        dx(:,1) = 0;
%         dx_ap = zeros(n,1);
%         dx_ap = dx(:,1) - diffX0;
        X(:,1) = X(:,1) + diffX0;
    end

end

