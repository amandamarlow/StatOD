function [GofXt, H] = GandH(t, X, stationNum, constants)
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

elevation = zeros(length(t),size(r_stn_E,2));
for i = 1:length(t) 
    r_sc_N = X(1:3, i);
    
    alpha = t(i)*omegaE + theta0;
    EN = Euler3(alpha);
    NE = EN';
   
    for m = 1:size(r_stn_E, 2)
       r_stn_N = NE*r_stn_E(:,m);
       v_stn_N = tilde(omegaE_N)*r_stn_N;       
       S_stn = [r_stn_N; v_stn_N];
       
       R = r_sc_N - r_stn_N;
       elevation(i,m) = asind(dot(R,r_stn_N)/norm(R)/norm(r_stn_N));
       
       if  elevation(i,m) >= 10
           [range, rangeRate, H_range, H_rangeRate] = Hcalcs(X(1:6, i), S_stn, NE);
           
           GofXt = [range; rangeRate];
           H = [H_range(1:n); H_rangeRate(1:n)];
       else
           GofXt = NaN(2,1);
           H = NaN(2,6);
       end
    end
end

% if measurementType == "range"
%     GofXt = range;
%     H = H_range;
% elseif measurementType == "rangeRate"
%     GofXt = rangeRate;
%     H = H_rangeRate;
% else
%     error("not a valid measurement type")
% end

end

