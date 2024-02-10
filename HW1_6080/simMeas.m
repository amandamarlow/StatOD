function [range_observations, rangeRate_observations, elevations_byStation, elevations_all] = simMeas(t, S_sc, r_stns_E, const)
%SIMMEAS Summary of this function goes here
%   Detailed explanation goes here

omegaE = const.omegaE;
omegaE_N = [0;0;omegaE];
theta0 = const.theta0;

% y = NaN(size(r_stns_E, 2), length(t));
n = 1;
elevations_all = zeros(length(t),size(r_stns_E,2));
for i = 1:length(t) 
    r_sc_N = S_sc(1:3, i);
    
    alpha = t(i)*omegaE + theta0;
    EN = Euler3(alpha);
    NE = EN';
   
    for m = 1:size(r_stns_E, 2)
       r_stn_N = NE*r_stns_E(:,m);
       v_stn_N = tilde(omegaE_N)*r_stn_N;       
       S_stn = [r_stn_N; v_stn_N];
       
       R = r_sc_N - r_stn_N;
%        elevation(i,m) = acosd(dot(R,r_stn_N)/norm(R)/norm(r_stn_N));
       elevations_all(i,m) = asind(dot(R,r_stn_N)/norm(R)/norm(r_stn_N));
       
       if  elevations_all(i,m) >= 10
%            [rng,rngRate] = Rng_RngRate(S_sc, S_stn);
%            y((m-1)*2+1:m*2, i) = [rng,rngRate];
           [range, rangeRate, H_range, H_rangeRate] = Hcalcs(S_sc(1:6, i), S_stn, NE);
           range_observations(n, :) = [t(i), m, range, H_range];
           rangeRate_observations(n, :) = [t(i), m, rangeRate, H_rangeRate];
           elevations_byStation(n,:) = [t(i), m, elevations_all(i,m)];
           n = n+1;
       end
    end
end

end

