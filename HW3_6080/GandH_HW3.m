function [GofXt, H] = GandH_HW3(t, X, stationNum, constants)
%SIMMEAS Summary of this function goes here
%   Detailed explanation goes here

    n = length(X);

    omegaE = constants.omegaE;
    omegaE_N = [0;0;omegaE];
    theta0 = constants.theta0;
    ae = constants.ae;

    % calculate station position
    stations_latlon = [
    -35.39833, 148.981944; 
    40.427222,  355.749444; 
    35.247164, 243.205
    ]*pi/180; % [rad] stations correspond to rows

    r_stn_E = latlon2ECEF(ae, stations_latlon(stationNum,1), stations_latlon(stationNum,2));
    
    alpha = t*omegaE + theta0;
    NE = [cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];
   [range, rangeRate, H_range, H_rangeRate] = Hcalcs_HW3(X(1:6), r_stn_E, NE, omegaE_N);

   GofXt = [range; rangeRate];

   H = [H_range(1:6); H_rangeRate(1:6)];
end

