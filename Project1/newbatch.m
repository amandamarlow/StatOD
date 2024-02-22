function [X, dx0, P, y, alpha, iterations, RMSresidual] = newbatch(t, data, R, X0, dx0, P0, constants)
%Batch filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
    t = data(:,1);
    n = length(X0);
    X = zeros(n, length(t));
    P = zeros(n, n, length(t));
    STM = zeros(n, n, length(t));
    H = NaN(2, n, length(t));
    y = NaN(2, length(t)); % pre-fit residual
    alpha = NaN(2, length(t)); % post-fit residual
    
    X(:,1) = X0;
%     P(:,:,1) = P0;
%     STM(:,:,1) = eye(n);
    P0est = P0;
    
    diffx0(:,1) = ones(n,1);
    max_iterations = 8;
    for q = 1:max_iterations  
        %% Initialize Iteration
        % or calculate rms of prefit and postfit residuals every loop and
        if min(data(:,1)) == t(1)
            Y = data(1,3:4)';
            stationNum = data(1,2);
            [GofXt0, H(:,:,1)] = GandH(t(1), X(:,1), stationNum, constants); % STM at t0 is I
            A = P0\eye(n) + H(:,:,1).'*(R(:,:,1)\H(:,:,1)); % Information Matrix
            y(:,1) = Y - GofXt0; % pre-fit residual
            N = H(:,:,1).'*(R(:,:,1)\y(:,1));       
        else
            A = P0\eye(n); % Information Matrix
            N = P0\dx0(:,q);
        end
        [X, STM] = integrateTrajectorySTM(data(:,1), X(:,1), eye(n), constants);
        for i = 2:length(data(:,1))
            Y = data(i,3:4)';
            stationNum = data(i,2);
            [GofXt, Htilde] = GandH(t(i), X(:,i), stationNum, constants);
            H(:,:,i) = Htilde*STM(:,:,i);
            y(:,i) = Y - GofXt; % pre-fit residual
            A = A + H(:,:,i).'*(R(:,:,i)\H(:,:,i));
            N = N + H(:,:,i).'*(R(:,:,i)\y(:,i));
            alpha(:,i) = y(:,i) - H(:,:,i)*diffx0(:,end); % post-fit residual
        end
        RMSresidual(:,q) = sqrt(1/length(t)*sum([alpha(1,~isnan(alpha(1,:)));alpha(2,~isnan(alpha(2,:)))] .^2, 2));
        % has the process converged?
%         if all(abs(diffx0(:,q))./(diag(P0(:,:,1).^(1/2))) <= 0.01)
        if all(abs(diffx0(:,q))./(diag(P0est).^1/2) <= 0.1)
            % should probably converge in 3 for project
            iterations = q;
            break
        elseif q == max_iterations
            iterations = q;
            warning("Failed to converge in %d iterations", max_iterations);
        end
        % solve normal equations
        diffx0(:,q+1) = A\N;
        dx0(:,q+1) = dx0(:,q) - diffx0(:,q+1);
        P0est = A\eye(n);
        X(:,1) = X(:,1) + diffx0(:,q+1);
        
    end
    for i = 2:length(t)
        P(:,:,i) = STM(:,:,i)*P0est*STM(:,:,i).';
    end
end

