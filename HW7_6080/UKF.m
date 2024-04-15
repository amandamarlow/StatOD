function [x, P, preFit] = UKF(t, data, R, Qc, x0, P0, alpha, beta, constants)
% Unscented Kalman Filter

    % L = length(x0)+length(Qc)+length(R);
    % kapa = 3-L;
    % Preallocate
    n = length(x0);
    x = zeros(n, length(t));
    P = zeros(n, n, length(t));
    % x_aug = zeros(L, length(t));
    % P_aug = zeros(L, L, length(t));
    preFit = zeros(2, length(t)); % pre-fit residual
    % postFit = zeros(2, length(t)); % post-fit residual
    %% Initialize
    x(:,1) = x0;
    P(:,:,1) = P0;
    % x_aug(:,1) = [x0; zeros(length(Qc),1); zeros(length(R),1)];
    % P_aug(:,:,1) = blkdiag(P0, Qc, R);
    measIdx = 1;
    % X = zeros(2*L+1, length(t))
    
    %% attempt at the augmented state vector
    % lambda = alpha^2*(L+kapa)-L;
    % 
    % X = zeros(L,2*L+1);
    % 
    % S = chol(P_aug(:,:,1));
    % X(1,:) = x_aug(:,1);
    % for j = 1:L
    %     X(j+1,:) = x_aug(:,1) + sqrt(L + lambda)*S(j,:)';
    %     X(L+j+1,:) = x_aug(:,1) - sqrt(L + lambda)*S(j,:)';
    % end
    % 
    % Wm = zeros(1,2*L+1);
    % Wm(2:2*L+1,:) = 1/(2*(L+lambda));
    % Wc = Wm;
    % Wm(1,:) = lambda/(L+lambda);
    % Wc(1,:) = lambda/(L+lambda) + (1-alpha^2+beta);
    
    kapa = 3-n;
    lambda = alpha^2*(n+kapa)-n;
    
    X = zeros(n,2*n+1, length(t));

    % S = chol(P(:,:,1));
    % X(:,1,1) = x(:,1);
    % for j = 1:n
    %     X(:,j+1, 1) = x(:,1) + sqrt(n + lambda)*S(j,:)';
    %     X(:,n+j+1, 1) = x(:,1) - sqrt(n + lambda)*S(j,:)';
    % end
    S = chol(P0);
    X(:,1,1) = x0;
    for j = 1:n
        X(:,j+1, 1) = x0 + sqrt(n + lambda)*S(:,j);
        X(:,n+j+1, 1) = x0 - sqrt(n + lambda)*S(:,j);
    end
    
    Wm = zeros(2*n+1, 1);
    Wm(2:end) = 1/(2*(n+lambda))*ones(2*n,1);
    Wc = Wm;
    Wm(1) = lambda/(n+lambda);
    Wc(1) = lambda/(n+lambda) + (1-alpha^2+beta);

    for k = 2:length(t)
        %% Time update
        X_ap = zeros(n,n*2+1);
        for i = 1:2*n+1
            X_ap(:,i) = F_UKF(t(k-1), t(k), X(:,i,k-1), constants);
        end
        % X(:,:,k) = zeros(n,n*2+1);
        % for i = 1:2*n+1
        %     X(:,i,k) = F_UKF(t(k-1), t(k), X(:,i,k-1), constants);
        % end

        % SNC
        dt = t(k) - t(k-1);
        if t(k) - data(measIdx,1) > 1000
            Qd = zeros(n);
        else
            Qd = [dt^3/3*Qc, dt^2/2*Qc;
                dt^2/2*Qc, dt*Qc];
        end
        
        x_ap = zeros(n,1);
        Psum = zeros(n,n);
        for i = 1:2*n+1
            x_ap = x_ap + Wm(i)*X_ap(:,i);
            % Psum = Psum + Wc(i)*(X_ap(:,i)-x_ap)*(X_ap(:,i)-x_ap)'; 
        end
        for i = 1:2*n+1
            % x_ap = x_ap + Wm(i)*X_ap(:,i);
            Psum = Psum + Wc(i)*(X_ap(:,i)-x_ap)*(X_ap(:,i)-x_ap)'; 
        end
        % x_ap = zeros(n,1);
        % Psum = zeros(n,n);
        % for i = 1:2*n+1
        %     x_ap = x_ap + Wm(i)*X(:,i);
        %     % Psum = Psum + Wc(i)*(X_ap(:,i)-x_ap)*(X_ap(:,i)-x_ap)'; 
        % end
        % for i = 1:2*n+1
        %     % x_ap = x_ap + Wm(i)*X_ap(:,i);
        %     Psum = Psum + Wc(i)*(X(:,i,k)-x_ap)*(X(:,i,k)-x_ap)'; 
        % end
        P_ap = Qd + Psum;
        
        S = chol(P_ap);
        X(:,1,k) = x_ap;
        for j = 1:n
            X(:,j+1, k) = x_ap + sqrt(n + lambda)*S(:,j);
            X(:,n+j+1, k) = x_ap - sqrt(n + lambda)*S(:,j);
        end         

        %% Measurement Update
        Y_at_t = find(t(k) == data(:,1));
        if isempty(Y_at_t)
            x(:,k) = x_ap;
            P(:,:,k) = P_ap;
        else
            % Resample
            % S = chol(P_ap);
            % X(:,1,k) = x_ap;
            % for j = 1:n
            %     X(:,j+1, k) = x_ap + sqrt(n + lambda)*S(:,j);
            %     X(:,n+j+1, k) = x_ap - sqrt(n + lambda)*S(:,j);
            % end
            % Process all measurements
            for measIdx = Y_at_t
                y = data(j,3:4)';
                stationNum = data(j,2);
                Y = zeros(2,2*n+1);
                y_ap = zeros(2,1);
                for i = 1:2*n+1
                    [Y(:,i), ~] = GandH_HW3(t(k), X(:,i,k), stationNum, constants);
                    y_ap = y_ap + Wm(i)*Y(:,i);
                end
                Pyy = zeros(2);
                Pxy = zeros(6,2);
                for i = 1:2*n+1
                    Pyy = Pyy + Wc(i)*(Y(:,i)-y_ap)*(Y(:,i)-y_ap)';
                    Pxy = Pxy + Wc(i)*(X(:,i,k)-x_ap)*(Y(:,i)-y_ap)';
                end
                Pyy = Pyy + R;
                preFit(:,k) = y-y_ap; % pre-fit residual
                K = Pxy/Pyy; % kalman gain
                x(:,k) = x_ap + K*preFit(:,k); % CHECK SIGNS
                P(:,:,k) = P_ap - K*Pyy*K';
                % postFit(:,k) = preFit(:,k) - H*dx; % post-fit residual 
            end
        end
    end

end

