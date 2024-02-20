function [GofXt, H] = GandH(t, X, stationNum, constants)
%SIMMEAS Summary of this function goes here
%   Detailed explanation goes here

    if stationNum == 101
        stationNum = 1;
    elseif stationNum == 337
        stationNum = 2;
    elseif stationNum == 394
        stationNum = 3;
    end

    n = length(X);

    omegaE = constants.omegaE;
    omegaE_N = [0;0;omegaE];
    theta0 = constants.theta0;
    ae = constants.ae;

    r_stn_E = X(10+3*(stationNum-1):9+3*(stationNum));

    r_sc_N = X(1:3);
    
    alpha = t*omegaE + theta0;
    EN = Euler3(alpha);
    NE = EN';
   

   r_stn_N = NE*r_stn_E;
   v_stn_N = tilde(omegaE_N)*r_stn_N;       
   S_stn = [r_stn_N; v_stn_N];

   R = r_sc_N - r_stn_N;
   elevation = asind(dot(R,r_stn_N)/norm(R)/norm(r_stn_N));

%    if  elevation >= 10
       [range, rangeRate, H_range, H_rangeRate] = Hcalcs(X, S_stn, NE);

       GofXt = [range; rangeRate];
       
       H = zeros(2, 18);
       H(:,1:6) = [H_range(1:6); H_rangeRate(1:6)];
       H(:,10+3*(stationNum-1):9+3*(stationNum)) = [H_range(7:9); H_rangeRate(7:9)];
%    else
%        GofXt = NaN(2,1);
%        H = NaN(2,n);
%    end

end

