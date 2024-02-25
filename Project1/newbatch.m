function [X, dx0, P, y, alpha, iterations, RMSresidual_post] = newbatch(data, R, X0, dx0, P0, measurement_type, constants)
%Batch filter
%   data -> [time after epoch [s], station number,
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
    if measurement_type == "range"
        H = zeros(1, n, length(t));
        y = zeros(1, length(t)); % pre-fit residual
        alpha = zeros(1, length(t)); % post-fit residual
        R = R(1,1);
    elseif measurement_type == "rangeRate"
        H = zeros(1, n, length(t));
        y = zeros(1, length(t)); % pre-fit residual
        alpha = zeros(1, length(t)); % post-fit residual
        R = R(2,2);
    else
        H = zeros(2, n, length(t));
        y = zeros(2, length(t)); % pre-fit residual
        alpha = zeros(2, length(t)); % post-fit residual
    end
    X(:,1) = X0;
    P0est = P0;

    diffx0 = zeros(n,1);
    max_iterations = 10;
    for q = 1:max_iterations
        %% Initialize Iteration
        if min(data(:,1)) == t(1)
            Y = data(1,3:4)';
            stationNum = data(1,2);
            [GofXt0, Hfull] = GandH(t(1), X(:,1), stationNum, constants); % STM at t0 is I
            if measurement_type == "range"
                Y = data(1,3);
                GofXt0 = GofXt0(1);
                H(:,:,1) = Hfull(1,:,1);
            elseif measurement_type == "rangeRate"
                Y = data(1,4);
                GofXt0 = GofXt0(2);
                H(:,:,1) = Hfull(2,:,1);
            elseif measurement_type == "both"
                H(:,:,1) = Hfull;
            end
            A = P0\eye(n) + H(:,:,1).'*(R\H(:,:,1)); % Information Matrix
            y(:,1) = Y - GofXt0; % pre-fit residual
            N = H(:,:,1)'*(R\y(:,1));   
            alpha(:,1) = y(:,1) - H(:,:,1)*diffx0; % post-fit residual
        else
            A = P0\eye(n); % Information Matrix
            N = P0\dx0;
        end
        [X, STM] = integrateTrajectorySTM(data(:,1), X(:,1), eye(n), constants);
        for i = 2:length(data(:,1))
            Y = data(i,3:4)';
            stationNum = data(i,2);
            [GofXt, Htilde] = GandH(t(i), X(:,i), stationNum, constants);
            if measurement_type == "range"
                Y = data(i,3);
                GofXt = GofXt(1);
                Htilde = Htilde(1,:);
            elseif measurement_type == "rangeRate"
                Y = data(i,4);
                GofXt = GofXt(2);
                Htilde = Htilde(2,:);
            end
            H(:,:,i) = Htilde*STM(:,:,i);
            y(:,i) = Y - GofXt; % pre-fit residual
            A = A + H(:,:,i)'*(R\H(:,:,i));
            N = N + H(:,:,i)'*(R\y(:,i));
            alpha(:,i) = y(:,i) - H(:,:,i)*diffx0; % post-fit residual
        end
%         RMSresidual_qminus = RMSresidual;
        RMSresidual_pre = sqrt(1/length(t)*sum(y.^2, 2));
        RMSresidual_post = sqrt(1/length(t)*sum(alpha.^2, 2));
%         if all(abs(RMSresidual-RMSresidual_qminus) < 0.01*RMSresidual)
        if all(abs(RMSresidual_post-RMSresidual_pre) < 0.01*RMSresidual_post) && q>1
            iterations = q;
            break
        elseif q == max_iterations && q>1
            iterations = q;
            warning("Failed to converge in %d iterations", max_iterations);
            break
        end
        % solve normal equations
        diffx0 = A\N;
        dx0 = dx0 - diffx0;
        P0est = A\eye(n);
        X(:,1) = X(:,1) + diffx0;
        
    end
    for i = 2:length(t)
        P(:,:,i) = STM(:,:,i)*P0est*STM(:,:,i)';
    end
end

