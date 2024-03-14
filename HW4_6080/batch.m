function [X, dx0, P, y, alpha, iterations, RMSresidual] = batch(t, Ymat, R, X0, dx0, P0, constants)
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
    P = zeros(n, n, length(t));
    STM = zeros(n, n, length(t));
    H = NaN(2, n, length(t));
    y = NaN(2, length(t)); % pre-fit residual
    alpha = NaN(2, length(t)); % post-fit residual
    
    X(:,1) = X0;
    P(:,:,1) = P0;
    STM(:,:,1) = eye(n);
    
    diffx0(:,1) = ones(n,1);
    max_iterations = 8;
    for q = 1:max_iterations  
        % change this to finite number or current diff - prev diff over covariance < 0.01 sigma
        %% Initialize Iteration
%         STM = eye(n);
        % or calculate rms of prefit and postfit residuals every loop and
        if min(Ymat(:,1)) == t(1)
            Y = Ymat(1,3:4)';
            stationNum = Ymat(1,2);
            [GofXt0, H(:,:,1)] = GandH_HW3(t(1), X(:,1), stationNum, constants); % STM at t0 is I
            A = P0(:,:,1)^-1 + H(:,:,1).'*(R(:,:,1)^-1)*H(:,:,1); % Information Matrix
            y(:,1) = Y - GofXt0; % pre-fit residual
            N = H(:,:,1).'*(R(:,:,1)^-1)*y(:,1);       
        else
            A = P0(:,:,1)^-1; % Information Matrix
            N = P0(:,:,1)^-1*dx0(:,q);
        end
        % compare and stop if small
        for i = 2:length(t)
            %% Integrate reference trajectory and STM from t(i-1) to t(i)
            [X(:,i), STM(:,:,i)] = integrateTrajectorySTM_HW3(t(i-1), t(i), X(:,i-1), STM(:,:,i-1), constants);
            %% Time update
%             dx_ap = STM*dx(:,i-1); % a priori state deviation estimate
%             P_ap = STM*P(:,:,i-1)*STM';
            %% Compute observation deviation, observation-state matrix, kalman gain matrix
            Y_at_t = find(t(i) == Ymat(:,1));
            if isempty(Y_at_t)
%                 y(:,i) = NaN(2,1); % pre-fit residual
%                 alpha(:,i) = NaN(2,1); % post-fit residual
            else
                for j = Y_at_t
                    Y = Ymat(j,3:4)';
                    stationNum = Ymat(j,2);
                    [GofXt, Htilde] = GandH_HW3(t(i), X(:,i), stationNum, constants);
                    H(:,:,i) = Htilde*STM(:,:,i);
                    if isnan(GofXt(1))
%                         y(:,i) = NaN(2,1); % pre-fit residual
%                         alpha(:,i) = NaN(2,1); % post-fit residual 
                    else
                        y(:,i) = Y - GofXt; % pre-fit residual
                        A = A + H(:,:,i).'*(R(:,:,j)^-1)*H(:,:,i);
                        N = N + H(:,:,i).'*(R(:,:,j)^-1)*y(:,i);
                    end
                end  
            end
            alpha(:,i) = y(:,i) - H(:,:,i)*diffx0(:,end); % post-fit residual
        end
        RMSresidual(:,q) = sqrt(1/length(t)*sum([alpha(1,~isnan(alpha(1,:)));alpha(2,~isnan(alpha(2,:)))] .^2, 2));
        % has the process converged?
%         if all(abs(diffx0(:,q))./(diag(P0(:,:,1).^(1/2))) <= 0.01)
        if all(abs(diffx0(:,q))./(diag(P0(:,:,q)).^1/2) <= 0.1)
            % should probably converge in 3 for project
            iterations = q;
            break
        elseif q == max_iterations
            iterations = q;
            warning("Failed to converge in %d iterations", max_iterations);
        end
        % solve normal equations
        diffx0(:,q+1) = inv(A)*N;
        dx0(:,q+1) = dx0(:,q) - diffx0(:,q+1);
        P0(:,:,q+1) = A^-1;
        X(:,1) = X(:,1) + diffx0(:,q+1);
        
    end
    for i = 2:length(t)
%         [X(:,i), STM(:,:,i)] = integrateTrajectorySTM(t(i-1), t(i), X(:,i-1), STM(:,:,i-1), constants);
        P(:,:,i) = STM(:,:,i)*P0(:,:,end)*STM(:,:,i).';
%         P(:,:,i) = STM(:,:,i)*P0(:,:,1)*STM(:,:,i).';
%         alpha(:,i) = y(:,i) - H(:,:,i)*dx0(:,end); % post-fit residual
%         alpha(:,i) = y(:,i) - H(:,:,i)*diffx0(:,end); % post-fit residual
    end
end
