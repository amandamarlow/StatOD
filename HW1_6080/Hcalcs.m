function [range, rangeRate, H_range, H_rangeRate] = Hcalcs(SC_state, station_state, ECEF2ECI)
%Hcalcs H matrix of range / range rate measurement partials with 
%respect to spacecraft state and station location
%
%   SC_state and station_state are the inertial frame position and velocity
%   states of the spacecraft and station respectively. 
%
%   Outputs the H matrix for range and range rate measurements with respect to
%   the spacecraft state and station location

R = SC_state(1:3);
V = SC_state(4:6);
Rs = station_state(1:3);
Vs = station_state(4:6);

range = norm(R-Rs);
rangeRate = dot((R-Rs),(V-Vs))/range;

% spacecraft state partial
rangePartialR = (R-Rs)'/range;
rangePartialV = zeros(1,3);
rangeRatePartialR = (V-Vs)'/range - rangeRate/(range^2)*(R-Rs)';
rangeRatePartialV = (R-Rs)'/range;

% station location partials
rangePartialRs = (Rs-R)'/range*ECEF2ECI;
rangeRatePartialRs = (-(V-Vs)'/range + rangeRate/(range^2)*(R-Rs)')*ECEF2ECI;

H_range = [rangePartialR, rangePartialV, rangePartialRs];
H_rangeRate = [rangeRatePartialR, rangeRatePartialV, rangeRatePartialRs];
end

