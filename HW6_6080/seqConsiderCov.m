function [X, Xref, dx, P, Pxx, Pxc, y, alpha] = seqConsiderCov(t, data, meas_cov, X0, dx0, P0, Pxx0, dc, Pcc, constants)
%Classic/Linearized Kalman Filter
%   Y = Measurement Matrix -> [time after epoch [s], station number,
%   range, rangeRate]
%   X0 = initial full state
%   dx0 = initial a priori estimated state deviation
%   P0 = error covatiance associated with dx0

    % Preallocate
%     t = Ymat(:,1); % time after epoch [s] corresponds to measurement vector
    n = length(X0);
    q = 1;
    Xref = zeros(n, length(t));
    dx = zeros(n, length(t));
    P = zeros(n, n, length(t));
    Pxx = zeros(n, n, length(t));
    Pxc = zeros(n, q, length(t));
    % Pcc = zeros(q, q, length(t));
    S = zeros(n, q, length(t));
    y = NaN(2, length(t)); % pre-fit residual
    alpha = NaN(2, length(t)); % post-fit residual
    j = 1;
    %% Initialize
    Xref(:,1) = X0;
    % dx(:,1) = dx0;
    P(:,:,1) = P0;
    Pxx(:,:,1) = Pxx0;
    S0_ap = zeros(n,q);
    Hc = zeros(2,q);

    %% Integrate reference trajectory and STM from t(i-1) to t(i)
    [Xref, PSI] = integrateTrajectoryPSI_consider(t, Xref(:,1), eye(n), zeros(n,q), constants);
    
    if data(1,1) == 0
        stationNum = data(1,2);
        [GofXt, Hx] = GandH_HW3(t(1), Xref(:,1), stationNum, constants);
        Y = data(1,3:4)';
        y(:,1) = Y - GofXt; % pre-fit residual
        [P(:,:,1), S(:,:,1), Pxx(:,:,1), Pxc, dx(:,1), dxc] = measUpdate_consider(Pcc, dc, P0, S0_ap, dx0, Hx, Hc, y(:,1), meas_cov);
        
        alpha(:,1) = y(:,1) - Hx*dx(:,1); % post-fit residual 
    end

    for k = 2:length(t)
        %% Time update
        PSI_tk_tkminus = PSI(:,:,k)/PSI(:,:,k-1);
        STM = PSI_tk_tkminus(1:n,1:n);
        theta = PSI_tk_tkminus(1:n,end-q+1:end);
        % [Xref(:,k), STM, theta] = integrateTrajectorySTM_consider([t(k-1), t(k)], Xref(:,k-1), eye(n), zeros(n,q), constants);
        [P_ap, S_ap, dx_ap, dxc_ap, Pxx_ap, Pxc_ap] = timeUpdate_consider(Pcc, dc, P(:,:,k-1), dx(:,k-1), S(:,:,k-1), STM, theta);
        %% Compute observation deviation, observation-state matrix, kalman gain matrix
        Y_at_t = find(t(k) == data(:,1));
        if isempty(Y_at_t)
            dx(:,k) = dx_ap;
            P(:,:,k) = P_ap;
            Pxx(:,:,k) = Pxx_ap;
            Pxc(:,:,k) = Pxc_ap;
            S(:,:,k) = S_ap;
        else
            for j = Y_at_t
                Y = data(j,3:4)';
                stationNum = data(j,2);
                [GofXt, Hx] = GandH_HW3(t(k), Xref(:,k), stationNum, constants);
                y(:,k) = Y - GofXt; % pre-fit residual
                [P(:,:,k), S(:,:,k), Pxx(:,:,k), Pxc(:,:,k), dx(:,k), dxc] = measUpdate_consider(Pcc, dc, P_ap, S_ap, dx_ap, Hx, Hc, y(:,k), meas_cov);
                alpha(:,k) = y(:,k) - Hx*dx_ap; % post-fit residual 
            end  
        end
    end

    X = Xref + dx;

end
