function [X, P, u, Xref, dx, ytilde, e] = SRIF_procNoise(t, data, X0, dx0, P0, meas_cov, Qc, constants)
%SRIF_BASIC Summary of this function goes here
%   Detailed explanation goes here

n = length(X0);
tsteps = length(t);
% Preallocate
ytilde = NaN(2, tsteps);
e = NaN(2, tsteps);
dx = zeros(n,tsteps);
u = zeros(n,tsteps-1);
P = zeros(n,n,tsteps);

%% Initialize
info = P0\eye(n);
R_ap = chol(info);
b_ap = R_ap*dx0;

%% Integrate ref trajectory
[Xref, STM_t0] = integrateTrajectorySTM_HW3(t, X0, eye(n), constants);

[dx(:,1), R] = measUpdateSRIF(t(1), data(1,:), X0, R_ap, b_ap, meas_cov, constants);
P(:,:,1) = R\inv(R');
for k = 2:tsteps
    %% SNC
    dt = t(k) - t(k-1);
    Qd = [dt^3/3*Qc, dt^2/2*Qc;
        dt^2/2*Qc, dt*Qc];
    [gamma,D,~] = svd(Qd);
    Ru = chol(inv(Qd));
    %% Time Update
    STM = STM_t0(:,:,k)/STM_t0(:,:,k-1);
    Rtilde = R/STM;

    % A = [Ru, zeros(n), bu_ap(:,k-1); -Rtilde*gamma, Rtilde, b(:,k-1)];
    A = [Ru, zeros(n), bu_ap; -Rtilde*gamma, Rtilde, b];
    [~,TA] = qr(A);
    Ru_ap = TA(1:n, 1:n);
    Rux_ap = TA(n+1:2*n,1:n);
    R_ap = TA(n+1:2*n, n+1:2*n);
    butilde = TA(1:n, end);
    b_ap = TA(n+1:end, end);
    
    x_ap = R_ap\b_ap;
    u_ap = Ru_ap\(bu-Ru*x_ap);
    
    %% Measurement Update
    Y_at_t = find(t(k) == data(:,1));
    if isempty(Y_at_t)
        dx(:,k) = x_ap;
        R = R_ap;
    else
        for j = Y_at_t
            [dx(:,k), R] = measUpdateSRIF(t(k), data(j,:), Xref(:,k), R_ap, b_ap, meas_cov, constants);
            R_ap = R;
            b_ap = R*dx(:,k);

        end  
    end
    % P(:,:,k) = (R'*R)\eye(n);
    u(:,k-1) = Ru_ap\(butilde-Ru_ap*dx(:,k));
    P(:,:,k) = R\inv(R');
end

X = Xref + dx;

end
