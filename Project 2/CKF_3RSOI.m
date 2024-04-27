function [t, X, Xref, dx, P, y, alpha] = CKF_3RSOI(t, data, R, Qc, Xref0, dx0, P0, iterations, smoothed, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    

%     % Preallocate
% %     t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(Xref0);
%     Xref = zeros(n, length(t));
%     dx = zeros(n, length(t));
%     P = zeros(n, n, length(t));
%     P_ap = zeros(n, n, length(t));
%     STM = zeros(n, n, length(t));
%     y = NaN(2, length(t)); % pre-fit residual
%     alpha = NaN(2, length(t)); % post-fit residual
%     j = 1;
%     %% Initialize
%     Xref(:,1) = Xref0;
%     % dx(:,1) = dx0;
%     P(:,:,1) = P0;
%     P_ap(:,:,1) = P0;
%     STM(:,:,1) = eye(n);
%     RMSresidual_pre = zeros(2,1);
    
    for q = 1:iterations
      
        %% Integrate reference trajectory and STM from t(i-1) to t(i)
        STM(:,:,1) = eye(n);
        [t, Xref, STM_t0] = integrateTo3RSOI(t, Xref0, STM(:,:,1), constants);

        % Preallocate
        dx = zeros(n, length(t));
        P = zeros(n, n, length(t));
        P_ap = zeros(n, n, length(t));
        % STM = zeros(n, n, length(t));
        y = NaN(2, length(t)); % pre-fit residual
        alpha = NaN(2, length(t)); % post-fit residual
        j = 1;
        %% Initialize
        Xref(:,1) = Xref0;
        % dx(:,1) = dx0;
        P(:,:,1) = P0;
        P_ap(:,:,1) = P0;
        STM(:,:,1) = eye(n);
        RMSresidual_pre = zeros(2,1);
        dx(:,1) = dx0;

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
            STM(:,:,i) = STM_t0(:,:,i)/STM_t0(:,:,i-1);
            dx_ap = STM(:,:,i)*dx(:,i-1); % a priori state deviation estimate
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
            P_ap(:,:,i) = STM(:,:,i)*P(:,:,i-1)*STM(:,:,i)' + Qd;
            P_ap(:,:,i) = 1/2*(P_ap(:,:,i) + P_ap(:,:,i)');
            %% Compute observation deviation, observation-state matrix, kalman gain matrix
            Y_at_t = find(t(i) == data(:,1));
            if isempty(Y_at_t)
                dx(:,i) = dx_ap;
                P(:,:,i) = P_ap(:,:,i);
            else
                for j = Y_at_t
                    for stnNum = 1:3
                        Y = data(j,[stnNum+1, stnNum+4])';
                        if ~isnan(Y(1))
                            [GofXt, H] = GandH_proj2(t(i), Xref(:,i), stnNum, constants);
                            y(:,i) = Y - GofXt; % pre-fit residual
                            K = P_ap(:,:,i)*H'/(H*P_ap(:,:,i)*H' + R); % kalman gain
                            dx(:,i) = dx_ap + K*(y(:,i) - H*dx_ap);
                            P(:,:,i) = (eye(n) - K*H)*P_ap(:,:,i)*(eye(n) - K*H)' + K*R*K';
                            alpha(:,i) = y(:,i) - H*dx(:,i); % post-fit residual
                        end
                    end
                end  
            end
        end
        diffx0 = STM_t0(:,:,end)\dx(:,end);
        % if all(dx(:,end) < diag(P(:,:,end)).^(1/2))
        RMSresidual_qminus = RMSresidual_pre;
        RMSresidual_pre = sqrt(1/length(t)*sum(y.^2, 2));
        if all(abs(RMSresidual_pre-RMSresidual_qminus) < 0.01*RMSresidual_pre) && q>1
            fprintf("converged in %d iterations \n", q)
            break
        elseif q == iterations && iterations > 1
            fprintf("failed to converge in %d iterations \n", iterations)
            break
        end
        Xref0 = Xref0+diffx0;
        dx0 = zeros(n,1);
    end

    if smoothed
        [dx, X, P] = smoother(t, Xref, dx, P, P_ap, STM);
    else
        X = Xref + dx;
    end
    
end
