function [range, rangeRate, H_range, H_rangeRate] = Hcalcs_proj2(X_sc, X_station, ECEF2ECI, omegaE_N)
%Hcalcs H matrix of range / range rate measurement partials with 
%respect to spacecraft state and station location
%   
%   X_station is the ECEF position of the station
%
%   X_sc is the inertial frame position and velocity
%   states of the spacecraft. +
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
rangeRate = r'*v/range;

% spacecraft state partial
% rangePartialR = (R-Rs)'/range;
rangePartialR = r'/range*eye(3);
rangePartialV = zeros(1,3);
% rangeRatePartialR = (V-Vs)'/range - rangeRate/(range^2)*(R-Rs)';
rangeRatePartialR = v'/range*(eye(3)-r/range*(r'*eye(3)/range));
rangeRatePartialV = r'/range*eye(3);

% station location partials
% rangePartialRs = -(R-Rs)'/range*ECEF2ECI;
rangePartialRs = -r'/range*ECEF2ECI*eye(3);
rangeRatePartialRs = -r'/range*tilde(omegaE_N)*ECEF2ECI*eye(3) + v'/range*(-ECEF2ECI*eye(3) - r/range*(r'/range*(-ECEF2ECI)*eye(3)));

H_range = [rangePartialR, rangePartialV, rangePartialRs];
H_rangeRate = [rangeRatePartialR, rangeRatePartialV, rangeRatePartialRs];
end

