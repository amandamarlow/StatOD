function [dx, R, ytilde] = measUpdateSRIF(t, data_k, Xref, R_ap, b_ap, meas_cov, constants)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(Xref);

Y = data_k(3:4)';
stationNum = data_k(2);
[GofXt, H] = GandH_HW3(t, Xref, stationNum, constants);
y = Y - GofXt; % pre-fit residual
%% Whitening
V = chol(meas_cov);
ytilde = V\y;
Htilde = V\H;

%% Householder's
for i = 1:2
    y_scalar = ytilde(i);
    A = [R_ap, b_ap; Htilde(i,:), y_scalar];
    [~, TA] = qr(A);
    R = TA(1:n,1:n);
    b = TA(1:6,end);
    % e(i,k) = TA(n+1:end, end);
    R_ap = R;
    b_ap = b;
end

%% backward substitution
for i = n:-1:1
    sum = 0;
    if i < n
        for j = i+1:n
            sum = sum + R(i,j)*dx(j);
        end
    end
    dx(i) = (b(i)-sum)/R(i,i);
end

end

