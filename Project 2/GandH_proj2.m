function [GofXt, H] = GandH_proj2(t, X, stnNum, const)
%SIMMEAS Summary of this function goes here
%   Detailed explanation goes here

    n = length(X);

    omegaE = const.omegaE;
    omegaE_N = [0;0;omegaE];
    theta0 = const.theta0;
    ae = const.ae;
    X_sc = X(1:6);
    
    alpha = t*omegaE + theta0;
    NE = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];

    if n > 7 
        r_stn_E = X(end-8+(stnNum-1)*(3):end-9+(stnNum)*(3));
        [range, rangeRate, H_range, H_rangeRate] = Hcalcs_proj2(X_sc, r_stn_E, NE, omegaE_N);
        H = [[H_range(1:6); H_rangeRate(1:6)], zeros(2, n-6)];
        H(:, end-8+(stnNum-1)*(3):end-9+(stnNum)*(3)) = [H_range(7:9);H_rangeRate(7:9)];
        % fprintf("using constants \n")
    else
        r_stn_E = latlon2ECEF(ae+const.DSS_alt(stnNum), const.DSS_latlon(stnNum,1), const.DSS_latlon(stnNum,2));
        [range, rangeRate, H_range, H_rangeRate] = Hcalcs_proj2(X_sc, r_stn_E, NE, omegaE_N);
        H = [[H_range(1:6); H_rangeRate(1:6)], zeros(2, n-6)];
    end

    GofXt = [range; rangeRate];
end

