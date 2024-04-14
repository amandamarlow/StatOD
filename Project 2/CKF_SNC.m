function [X, Xref, dx, P, y, alpha] = CKF_SNC(t, data, R, Qc, Xref0, dx0, P0, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
%     t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(Xref0);
    Xref = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    y = NaN(2, length(t)); % pre-fit residual
    alpha = NaN(2, length(t)); % post-fit residual
    j = 1;
    %% Initialize
    Xref(:,1) = Xref0;
    % dx(:,1) = dx0;
    P(:,:,1) = P0;
%     STM = eye(n);

    %% Integrate reference trajectory and STM from t(i-1) to t(i)
    [Xref, STM_t0] = integrateTrajectorySTM_proj2(t, Xref0, eye(n), constants);
    for stnNum = 1:3
        Y = data(1,[stnNum+1, stnNum+4])';
        if ~isnan(Y(1))
            [GofXt, H] = GandH_proj2(t(1), Xref0, stnNum, constants);
            y(:,1) = Y - GofXt; % pre-fit residual
            K = P0*H'/(H*P0*H' + R); % kalman gain
            dx(:,1) = dx0 + K*(y(:,1) - H*dx0);
            P(:,:,1) = (eye(n) - K*H)*P0*(eye(n) - K*H)' + K*R*K';
            alpha(:,1) = y(:,1) - H*dx(:,1); % post-fit residual 
        end
    end
    for i = 2:length(t)
        %% Time update
        STM = STM_t0(:,:,i)/STM_t0(:,:,i-1);
        dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
        %% SNC
        dt = t(i) - t(i-1);
        if t(i) - data(j,1) > 1000
            Qd = zeros(n);
        else
            Qd = zeros(n);
            Qd(1:6,1:6) = [dt^3/3*Qc, dt^2/2*Qc;
                dt^2/2*Qc, dt*Qc];
        end
        % Qd = [dt^3/3*Qc, dt^2/2*Qc;
        %     dt^2/2*Qc, dt*Qc];
        P_ap = STM*P(:,:,i-1)*STM' + Qd;
        P_ap = 1/2*(P_ap + P_ap');
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(i) == data(:,1));
        if isempty(Y_at_t)
            dx(:,i) = dx_ap;
            P(:,:,i) = P_ap;
        else
            for j = Y_at_t
                for stnNum = 1:3
                    Y = data(j,[stnNum+1, stnNum+4])';
                    if ~isnan(Y(1))
                        [GofXt, H] = GandH_proj2(t(i), Xref(:,i), stnNum, constants);
                        y(:,i) = Y - GofXt; % pre-fit residual
                        K = P_ap*H'/(H*P_ap*H' + R); % kalman gain
                        dx(:,i) = dx_ap + K*(y(:,i) - H*dx_ap);
                        P(:,:,i) = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*R*K';
                        alpha(:,i) = y(:,i) - H*dx(:,i); % post-fit residual
                    end
                end
            end  
        end
    end
    X = Xref + dx;
end

