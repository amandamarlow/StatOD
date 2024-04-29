function [X, P, u, Xref, dx, y, e] = SRIF_procNoise(t, data, X0, dx0, P0, meas_cov, Qc, constants)
%SRIF_BASIC Summary of this function goes here
%   Detailed explanation goes here

n = length(X0);
q = length(Qc);
tsteps = length(t);
% Preallocate
y = NaN(2, tsteps);
e = NaN(2, tsteps);
dx = zeros(n,tsteps);
u = zeros(q,tsteps-1);
P = zeros(n,n,tsteps);
bu_ap = zeros(q,tsteps);
Ru = zeros(q,q,tsteps);
% Ru_ap = zeros(n,n,tsteps);
% gamma = zeros(n,n,tsteps);

%% Initialize
info = P0\eye(n);
R_ap = chol(info);
b_ap = R_ap*dx0;

%% Integrate ref trajectory
[Xref, STM_t0] = integrateTrajectorySTM_proj2(t, X0, eye(n), constants);
%% 1st step
for stnNum = 1:3
        Y = data(1,[stnNum+1, stnNum+4])';
    if ~isnan(Y(1))
        [dx(:,1), R, y(:,1), e(:,1)] = measUpdateSRIF(t(1), [data(1,1),stnNum,Y'], Xref(:,1), R_ap, b_ap, meas_cov, constants);
    end
end
% [dx(:,1), R, y(:,1)] = measUpdateSRIF(t(1), data(1,:), X0, R_ap, b_ap, meas_cov, constants);
P(:,:,1) = R\inv(R');
b = R*dx(:,1);
for k = 2:tsteps
    %% SNC
    dt = t(k) - t(k-1);
    gamma = [dt^2/2*eye(3); dt*eye(3); zeros(n-6,3)];
    Ru(:,:,k) = chol(inv(Qc));
    %% Time Update
    STM = STM_t0(:,:,k)/STM_t0(:,:,k-1);
    Rtilde = R/STM;

    % A = [Ru, zeros(n), bu_ap(:,k-1); -Rtilde*gamma, Rtilde, b(:,k-1)];
    % A = [Ru(:,:,k-1), zeros(q,n), bu_ap(:,k-1); -Rtilde*gamma, Rtilde, b];
    A = [Ru(:,:,k), zeros(q,n), bu_ap(:,k-1); -Rtilde*gamma, Rtilde, b];
    [~,TA] = qr(A);
    % Ru_ap(:,:,k) = TA(1:n, 1:n);
    Ru_ap = TA(1:q, 1:q);
    % Rux_ap(:,:,k) = TA(n+1:2*n,1:n);
    % Rux_ap = TA(q+1:q+n,1:q);
    Rux_ap = TA(1:q, q+1:q+n);
    R_ap = TA(q+1:q+n, q+1:q+n);
    % butilde(:,k) = TA(1:n, end);
    butilde = TA(1:q, end);
    b_ap = TA(q+1:end, end);
    
    x_ap = R_ap\b_ap;
    u_ap = Ru_ap\(butilde-Rux_ap*x_ap);
    bu_ap(:,k-1) = Ru(:,:,k)*u_ap;
    % bu_ap(:,k) = Ru(:,:,k)*u_ap;
    % b_ap = R_ap*x_ap;
    %% Measurement Update
    Y_at_t = find(t(k) == data(:,1));
    if isempty(Y_at_t)
        dx(:,k) = R_ap\b_ap;
        R = R_ap;
        b = R*dx(:,k);
    else
        for j = Y_at_t
            for stnNum = 1:3
                    Y = data(j,[stnNum+1, stnNum+4])';
                if ~isnan(Y(1))
                    [dx(:,k), R, y(:,k), e(:,1)] = measUpdateSRIF(t(k), [data(j,1),stnNum,Y'], Xref(:,k), R_ap, b_ap, meas_cov, constants);
                    R_ap = R;
                    b = R*dx(:,k);
                    b_ap = b;
                end
            end
        end  
    end
    % P(:,:,k) = (R'*R)\eye(n);
    % u(:,k-1) = Ru_ap(:,:,k)\(butilde(:,k)-Rux_ap(:,:,k)*dx(:,k));
    u(:,k-1) = Ru_ap\(butilde-Rux_ap*dx(:,k));
    P(:,:,k) = R\inv(R');
    % P(:,:,k) = (R*R')\eye(n);
    % P(:,:,k) = (R'*R)\eye(n);
end

X = Xref + dx;

end

