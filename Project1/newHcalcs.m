function [range, rangeRate, H_range, H_rangeRate] = newHcalcs(X_sc, X_station, ECEF2ECI, omegaE_N)
%Hcalcs H matrix of range / range rate measurement partials with 
%respect to spacecraft state and station location
%
%   SC_state and station_state are the inertial frame position and velocity
%   states of the spacecraft and station respectively. 
%
%   Outputs the H matrix for range and range rate measurements with respect to
%   the spacecraft state and station location

R = X_sc(1:3);
V = X_sc(4:6);
Rs = X_station(1:3);
% Vs = X_station(4:6);

r = (R - ECEF2ECI*Rs);
v = (V - tilde(omegaE_N)*ECEF2ECI*Rs);

% range = norm(R-Rs);
range = sqrt(r'*r);
rangeRate = r'*v./range;

% spacecraft state partial
% rangePartialR = (R-Rs)'/range;
rangePartialR = r'*eye(3)./range;
rangePartialV = zeros(1,3);
% rangeRatePartialR = (V-Vs)'/range - rangeRate/(range^2)*(R-Rs)';
rangeRatePartialR = v'*(eye(3)-r*(r'*eye(3))./(range^2))./range;
rangeRatePartialV = r'/range*eye(3);

% station location partials
% rangePartialRs = -(R-Rs)'/range*ECEF2ECI;
rangePartialRs = -r'*ECEF2ECI*eye(3)./range;
rangeRatePartialRs = -r'*tilde(omegaE_N)*ECEF2ECI*eye(3)./range + v'*(-ECEF2ECI*eye(3) - r*(r'*(-ECEF2ECI)*eye(3)./range)./range)./range;

H_range = [rangePartialR, rangePartialV, rangePartialRs];
H_rangeRate = [rangeRatePartialR, rangeRatePartialV, rangeRatePartialRs];
end

