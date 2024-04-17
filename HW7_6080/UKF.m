function [x, P, preFit] = UKF(t, data, R, Qc, x0, P0, alpha, beta, constants)
% Unscented Kalman Filter
    % Preallocate
    n = length(x0);
    x = zeros(n, length(t));
    P = zeros(n, n, length(t));
    preFit = zeros(2, length(t)); % pre-fit residual
    % postFit = zeros(2, length(t)); % post-fit residual
    %% Initialize
    x(:,1) = x0;
    P(:,:,1) = P0;
    measIdx = 1;
  
    kapa = 3-n;
    lambda = alpha^2*(n+kapa)-n;
    
    X = zeros(n,2*n+1, length(t));

    S = chol(P0);
    X(:,1,1) = x0;
    for j = 1:n
        X(:,j+1, 1) = x0 + sqrt(n + lambda)*S(j,:)';
        X(:,n+j+1, 1) = x0 - sqrt(n + lambda)*S(j,:)';
    end
    
    Wm = zeros(2*n+1, 1);
    Wm(2:end) = 1/(2*(n+lambda))*ones(2*n,1);
    Wc = Wm;
    Wm(1) = lambda/(n+lambda);
    Wc(1) = lambda/(n+lambda) + (1-alpha^2+beta);

    for k = 2:length(t)
        % SNC
        dt = t(k) - t(k-1);
        if t(k) - data(measIdx,1) > 1000
            Qd = zeros(n);
        else
            Qd = [dt^3/3*Qc, dt^2/2*Qc;
                dt^2/2*Qc, dt*Qc];
        end
        
        for i = 1:(2*n+1)
            X(:,i,k) = F_UKF(t(k-1), t(k), X(:,i,k-1), constants);
        end

        x_ap = zeros(n,1);
        Psum = zeros(n,n);
        for i = 1:(2*n+1)
            x_ap = x_ap + Wm(i)*X(:,i,k);
        end
        for i = 1:(2*n+1)
            Psum = Psum + Wc(i)*(X(:,i,k)-x_ap)*(X(:,i,k)-x_ap)'; 
        end
        P_ap = Qd + Psum;
              

        %% Measurement Update
        Y_at_t = find(t(k) == data(:,1));
        if isempty(Y_at_t)
            x(:,k) = x_ap;
            P(:,:,k) = P_ap;
        else
            % Process all measurements
            for measIdx = Y_at_t
                y = data(measIdx,3:4)';
                stationNum = data(measIdx,2);
                Y = zeros(2,2*n+1);
                y_ap = zeros(2,1);
                for i = 1:(2*n+1)
                    [Y(:,i), ~] = GandH_HW3(t(k), X(:,i,k), stationNum, constants);
                    y_ap = y_ap + Wm(i)*Y(:,i);
                end
                Pyy = zeros(2);
                Pxy = zeros(6,2);
                for i = 1:(2*n+1)
                    Pyy = Pyy + Wc(i)*((Y(:,i)-y_ap)*(Y(:,i)-y_ap)');
                    Pxy = Pxy + Wc(i)*((X(:,i,k)-x_ap)*(Y(:,i)-y_ap)');
                end
                Pyy = Pyy + R;
                preFit(:,k) = y-y_ap; % pre-fit residual
                K = Pxy/Pyy; % kalman gain
                x(:,k) = x_ap + K*preFit(:,k);
                P(:,:,k) = P_ap - K*Pyy*K';
                % postFit(:,k) = preFit(:,k) - H*dx; % post-fit residual 
            end
        end
        S = chol(P(:,:,k));
        X(:,1,k) = x(:,k);
        for j = 1:n
            X(:,j+1, k) = x(:,k) + sqrt(n + lambda)*S(j,:)';
            X(:,n+j+1, k) = x(:,k) - sqrt(n + lambda)*S(j,:)';
        end   
    end

end

