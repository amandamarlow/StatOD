function [BdotVec, B_STM] = Bplane(r_N, v_N, mu)
%UNTITLED Summary of this function goes here
%   outputs BdotVec = [BdotR;BdotT]
    r = norm(r_N);
    v = norm(v_N);
    h_N = cross(r_N, v_N);
    h = norm(h_N);
    e_N = 1/mu*((v^2 - mu/r)*r_N - dot(r_N, v_N)*v_N);
    % e_N = cross(v_N,h_N)/mu - r_N/r;
    e = norm(e_N);
    P_hat = e_N/norm(e_N);
    energy = v^2/2-mu/r;
    a = -mu/2/energy;
    p = h^2/mu;
    b = abs(a)*sqrt(e^2 - 1);
    W_hat = h_N/h;
    Q_hat = cross(W_hat, P_hat);
    S_hat = 1/e*P_hat + sqrt(e^2-1)/e*Q_hat;

    % S_hat = v_N/v;
    N_hat = [0;0;1];
    % N_hat = [1;0;0];
    T_hat = cross(S_hat, N_hat)/norm(cross(S_hat, N_hat));
    R_hat = cross(S_hat, T_hat);

    B_N = b*cross(S_hat, W_hat);

    BdotR = dot(B_N, R_hat);
    BdotT = dot(B_N, T_hat);
    DCM_BN = [S_hat'; T_hat'; R_hat'];

    %partials
    vinf = sqrt(-mu/a);  
    % vinf = mu/h*sqrt(e^2-1);  
    alpha = vinf^2;

    v_part_r = zeros(3);
    v_part_v = eye(3);
    r_part_r = eye(3);
    r_part_v = zeros(3);

    h_part_r = zeros(3);
    h_part_v = zeros(3);
    % alpha_part_r = zeros(3);
    % alpha_part_v = zeros(3);
    for i = 1:3
        h_part_r(:,i) = cross(r_part_r(:,i), v_N) + cross(r_N, v_part_r(:,i));
        h_part_v(:,i) = cross(r_part_v(:,i), v_N) + cross(r_N, v_part_v(:,i));

    end

    alpha_part_r = 2*(v_N'*v_part_r + mu/(r^3)*r_N'*r_part_r);
    alpha_part_v = 2*(v_N'*v_part_v + mu/(r^3)*r_N'*r_part_v);
    a_part_r = -a/alpha*alpha_part_r;
    a_part_v = -a/alpha*alpha_part_v;
    vinf_part_r = alpha_part_r/2/vinf;
    vinf_part_v = alpha_part_v/2/vinf;

    BdotR_part_r = (T_hat'*h_part_r - BdotR*vinf_part_r)/vinf;
    BdotT_part_r = -(R_hat'*h_part_r + BdotT*vinf_part_r)/vinf;
    BdotR_part_v = (T_hat'*h_part_v - BdotR*vinf_part_v)/vinf;
    BdotT_part_v = -(R_hat'*h_part_v + BdotT*vinf_part_v)/vinf;

    % B_STM = [BdotR_part_r, BdotR_part_v; BdotT_part_r, BdotT_part_v];
    B_STM = [BdotT_part_r, BdotT_part_v; BdotR_part_r, BdotR_part_v];
    BdotVec = [BdotR; BdotT];
end
