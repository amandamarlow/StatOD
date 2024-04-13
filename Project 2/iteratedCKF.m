function [X, Xref, dx, P, y, alpha] = iteratedCKF(data, meas_cov, X0, dx0, P0, iterations, measurement_type, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
    t = data(:,1);
    n = length(X0);
    Xref = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    meas_cov = diag([meas_cov(1)*ones(1,3), meas_cov(2)*ones(1,3)]);
    % if measurement_type == "range"
    %     y = zeros(1, length(t)); % pre-fit residual
    %     alpha = zeros(1, length(t)); % post-fit residual
    %     R = R(1,1);
    % elseif measurement_type == "rangeRate"
    %     y = zeros(1, length(t)); % pre-fit residual
    %     alpha = zeros(1, length(t)); % post-fit residual
    %     R = R(2,2);
    % else
        y = zeros(6, length(t)); % pre-fit residual
        alpha = zeros(6, length(t)); % post-fit residual
    % end
    %% Initialize
    Xref(:,1) = X0;
    dx(:,1) = dx0;
%     dx_ap = dx0;
    P(:,:,1) = P0;
%     RMSresidual_post = zeros(2,1);
    RMSresidual_pre = zeros(6,1);
    
    for q = 1:iterations
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        [Xref, STM_t0] = integrateTrajectorySTM_proj2(t, Xref(:,1), eye(n), constants);
        stnsInView = ~isnan(data(2:4));
        [GofXt, H] = GandH_proj2(t(1), Xref(:,1), stnsInView, constants);
        % if measurement_type == "range"
        %     Y = data(1,3);
        %     GofXt = GofXt(1);
        %     H = H(1,:);
        % elseif measurement_type == "rangeRate"
        %     Y = data(1,4);
        %     GofXt = GofXt(2);
        %     H = H(2,:);
        % else
            Y = data(1,2:end)';
        % end
        y(:,1) = Y - GofXt; % pre-fit residual
%         K = P0*H'/(H*P0*H' + R); % kalman gain
%         dx(:,1) = dx_ap + K*(y(:,1) - H*dx_ap);
%         P(:,:,1) = (eye(n) - K*H)*P0*(eye(n) - K*H)' + K*R*K';
        alpha(:,1) = y(:,1) - H*dx(:,1); % post-fit residual 
        for i = 2:length(t)
            STM = STM_t0(:,:,i)/STM_t0(:,:,i-1);
            %% Time update
            dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
            P_ap = STM*P(:,:,i-1)*STM';

            Y = data(i,2:end)';
            stnsInView = ~isnan(Y(1:3));
            [GofXt, H] = GandH_proj2(t(i), Xref(:,i), stnsInView, constants);
            % if measurement_type == "range"
            %     Y = data(i,3);
            %     GofXt = GofXt(1);
            %     H = H(1,:);
            % elseif measurement_type == "rangeRate"
                % Y = data(i,4);
                % GofXt = GofXt(2);
                % H = H(2,:);
            % end
            y(:,i) = Y - GofXt; % pre-fit residual
            K = P_ap*H'/(H*P_ap*H' + meas_cov); % kalman gain
            dx_ap = dx_ap + K*(y(:,i) - H*dx_ap);
            P_ap = (eye(n) - K*H)*P_ap*(eye(n) - K*H)' + K*meas_cov*K';
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
        Xref(:,1) = Xref(:,1) + diffX0;
    end
    X = Xref + dx;
end

