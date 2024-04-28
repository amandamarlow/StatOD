function [BdotVec, B_STM] = Bplane(r_N, v_N, mu, varyShat)
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
    T_hat = cross(S_hat, N_hat)/norm(cross(S_hat, N_hat));
    % V_hat = [1;0;0];
    % T_hat = cross(V_hat, S_hat)/norm(cross(V_hat, S_hat));
    R_hat = cross(S_hat, T_hat);

    B_N = b*cross(S_hat, W_hat);

    BdotR = dot(B_N, R_hat);
    BdotT = dot(B_N, T_hat);
    DCM_BN = [S_hat'; T_hat'; R_hat'];

    %partials
    vinf = sqrt(-mu/a); 
    % vinf = v; 
    alpha = vinf^2;

    v_part_r = zeros(3);
    v_part_v = eye(3);
    r_part_r = eye(3);
    r_part_v = zeros(3);

    hN_part_r = zeros(3);
    hN_part_v = zeros(3);
    % alpha_part_r = zeros(3);
    % alpha_part_v = zeros(3);
    for i = 1:3
        hN_part_r(:,i) = cross(r_part_r(:,i), v_N) + cross(r_N, v_part_r(:,i));
        hN_part_v(:,i) = cross(r_part_v(:,i), v_N) + cross(r_N, v_part_v(:,i));
    end
    h_part_r = h_N'*hN_part_r/h;
    h_part_v = h_N'*hN_part_v/h;

    alpha_part_r = 2*(v_N'*v_part_r + mu/(r^3)*r_N'*r_part_r);
    alpha_part_v = 2*(v_N'*v_part_v + mu/(r^3)*r_N'*r_part_v);
    a_part_r = -a/alpha*alpha_part_r;
    a_part_v = -a/alpha*alpha_part_v;
    vinf_part_r = alpha_part_r/2/vinf;
    vinf_part_v = alpha_part_v/2/vinf;
    
    % Simpler Partials
    if ~varyShat
        BdotR_part_r = (T_hat'*hN_part_r - BdotR*vinf_part_r)/vinf;
        BdotT_part_r = -(R_hat'*hN_part_r + BdotT*vinf_part_r)/vinf;
        BdotR_part_v = (T_hat'*hN_part_v - BdotR*vinf_part_v)/vinf;
        BdotT_part_v = -(R_hat'*hN_part_v + BdotT*vinf_part_v)/vinf;
    elseif varyShat
        % Better Partials
    
        eN_part_r = zeros(3);
        eN_part_v = zeros(3);
        for i = 1:3
            eN_part_r(:,i) = (cross(v_part_r(:,i),h_N) + cross(v_N, hN_part_r(:,i))) - r_part_r(:,i)/r + dot(r_N,r_part_r(:,i))/r^3*r_N;
            eN_part_v(:,i) = (cross(v_part_v(:,i),h_N) + cross(v_N, hN_part_v(:,i))) - r_part_v(:,i)/r + dot(r_N,r_part_v(:,i))/r^3*r_N;
        end
        e_part_r = e_N'*eN_part_r/e;
        e_part_v = e_N'*eN_part_r/e;
    
        P_part_r = eN_part_r/e - e_part_r*e_N/e^2;
        P_part_v = eN_part_v/e - e_part_v*e_N/e^2;
    
        Q_part_r = zeros(3);
        Q_part_v = zeros(3);
        for i = 1:3
            Q_part_r(:,i) = (cross(hN_part_r(:,i),P_hat) + cross(h_N, P_part_r(:,i)) - h_part_r(:,i)*Q_hat)/h;
            Q_part_v(:,i) = (cross(hN_part_v(:,i),P_hat) + cross(h_N, P_part_v(:,i)) - h_part_v(:,i)*Q_hat)/h;
        end
    
        S_part_r = -e_part_r*P_hat/e + P_part_r/e + e_part_r*Q_hat/e^2/sqrt(e^2-1) + sqrt(e^2-1)*Q_part_r/e;
        S_part_v = -e_part_v*P_hat/e + P_part_v/e + e_part_v*Q_hat/e^2/sqrt(e^2-1) + sqrt(e^2-1)*Q_part_v/e;
        
        T_part_r = zeros(3,3);
        T_part_v = zeros(3,3);
        R_part_r = zeros(3,3);
        R_part_v = zeros(3,3);
        for i = 1:3
            % T_part_r(:,i) = cross(V_hat, S_part_r(:,i)) - T_hat'*cross(V_hat,S_part_r(:,i))*T_hat;
            % T_part_v(:,i) = cross(V_hat, S_part_v(:,i)) - T_hat'*cross(V_hat,S_part_v(:,i))*T_hat;
            T_part_r(:,i) = (cross(S_part_r(:,i), N_hat) - T_hat'*cross(S_part_r(:,i), N_hat)*T_hat)/norm(cross(S_hat,N_hat));
            T_part_v(:,i) = (cross(S_part_v(:,i), N_hat) - T_hat'*cross(S_part_v(:,i), N_hat)*T_hat)/norm(cross(S_hat,N_hat));
            
            R_part_r(:,i) = cross(S_part_r(:,i),T_hat) + cross(S_hat,T_part_r(:,i));
            R_part_v(:,i) = cross(S_part_v(:,i),T_hat) + cross(S_hat,T_part_v(:,i));
        end
    
        BdotR_part_r = (T_hat'*hN_part_r + h_N'*T_part_r - BdotR*vinf_part_r)/vinf;
        BdotT_part_r = -(R_hat'*hN_part_r + h_N'*R_part_r + BdotT*vinf_part_r)/vinf;
        BdotR_part_v = (T_hat'*hN_part_v + h_N'*T_part_v - BdotR*vinf_part_v)/vinf;
        BdotT_part_v = -(R_hat'*hN_part_v + h_N'*R_part_v + BdotT*vinf_part_v)/vinf;
    else
        error("Not a valit varyShat argument")
    end
    
    

    % B_STM = [BdotR_part_r, BdotR_part_v; BdotT_part_r, BdotT_part_v];
    B_STM = [BdotT_part_r, BdotT_part_v; BdotR_part_r, BdotR_part_v];
    % BdotVec = [BdotR; BdotT];
    BdotVec = [BdotT; BdotR];
end

