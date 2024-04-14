function [GofXt, H] = GandH_proj2(t, X, stnsInView, const)
%SIMMEAS Summary of this function goes here
%   Detailed explanation goes here

    n = length(X);

    omegaE = const.omegaE;
    omegaE_N = [0;0;omegaE];
    theta0 = const.theta0;
    ae = const.ae;
    X_sc = X(1:6);

    % % calculate station position
    % stations_latlon = [
    % -35.39833, 148.981944; 
    % 40.427222,  355.749444; 
    % 35.247164, 243.205
    % ]*pi/180; % [rad] stations correspond to rows
    
    alpha = t*omegaE + theta0;
    NE = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];

    GofXt = NaN(6,1);
    H = NaN(6,n);
    for stnNum = 1:length(stnsInView)
        if stnsInView(stnNum) == true
            r_stn_E = latlon2ECEF(const.DSS_alt(stnNum), const.DSS_latlon(stnNum,1), const.DSS_latlon(stnNum,2));
            [range, rangeRate, H_range, H_rangeRate] = Hcalcs_proj2(X_sc, r_stn_E, NE, omegaE_N);
            GofXt(stnNum) = range;
            GofXt(3+stnNum) = rangeRate;
            H(stnNum,:) = [H_range(1:6), zeros(1,n-6)];
            H(3+stnNum,:) = [H_rangeRate(1:6), zeros(1,n-6)];
        end
    end
    
   % [range, rangeRate, H_range, H_rangeRate] = Hcalcs_HW3(X(1:6), r_stn_E, NE, omegaE_N);
   % 
   % GofXt = [range; rangeRate];
   % 
   % H = [H_range(1:6); H_rangeRate(1:6)];
end

