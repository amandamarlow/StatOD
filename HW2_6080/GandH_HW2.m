function [GofXt, H] = GandH_HW2(t, X, stationNum, constants)
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

    r_sc_N = X(1:3);
    
    alpha = t*omegaE + theta0;
    EN = Euler3(alpha);
    NE = EN';
   

   r_stn_N = NE*r_stn_E;
   v_stn_N = tilde(omegaE_N)*r_stn_N;       
   S_stn = [r_stn_N; v_stn_N];

   R = r_sc_N - r_stn_N;
   elevation = asind(dot(R,r_stn_N)/norm(R)/norm(r_stn_N));

   if  elevation >= 10
       [range, rangeRate, H_range, H_rangeRate] = Hcalcs(X, S_stn, NE);

       GofXt = [range; rangeRate];
       H = [H_range(1:n); H_rangeRate(1:n)];
   else
       GofXt = NaN(2,1);
       H = NaN(2,6);
   end

end

