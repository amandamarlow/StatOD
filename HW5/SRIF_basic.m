function [X, P, Xref, dx, ytilde, e] = SRIF_basic(t, data, X0, dx0, P0, meas_cov, constants)
%SRIF_BASIC Summary of this function goes here
%   Detailed explanation goes here

n = length(X0);
tsteps = length(t);
% Preallocate
ytilde = NaN(2, tsteps);
e = NaN(2, tsteps);
dx = zeros(n,tsteps);
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
    %% Time Update
    STM = STM_t0(:,:,k)/STM_t0(:,:,k-1);
    R_ap = R/STM;
    [~, R_ap] = qr(R_ap);
    x_ap = STM*dx(:,k-1);
    b_ap = R_ap*x_ap;
    % b_ap = b; % WHY DOESN'T THIS WORK? - wrong sign
    
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
    P(:,:,k) = R\inv(R');
end

X = Xref + dx;

end

