function [range, rangeRate, H_range, H_rangeRate] = Hcalcs(SC_state, station_state, ECEF2ECI)
%RANGEPARTIALSCSTATE H matrix of range / range rate measurement partials with 
%respect to spacecraft state
%
%   SC_state and station_state are the inertial frame position and velocity
%   states of the spacecraft and station respectively. 
%
%   Outputs the H matrix for range and range rate measurements with respect to
%   the spacecraft state

R = SC_state(1:3);
V = SC_state(4:6);
Rs = station_state(1:3);
Vs = station_state(4:6);

x = R(1);
y = R(2);
z = R(3);
xdot = V(1);
ydot = V(2);
zdot = V(3);
xs = Rs(1);
ys = Rs(2);
zs = Rs(3);
xsdot = Vs(1);
ysdot = Vs(2);
zsdot = Vs(3);

range = norm(R-Rs);
rangeRate = dot((R-Rs),(V-Vs))/range;

% spacecraft state partial
rangePartialR = [x-xs, y-ys, z-zs]/range;
rangePartialV = zeros(1,3);
rangeRatePartialR = [xdot-xsdot, ydot-ysdot, zdot-zsdot]/range - rangeRate/(range^2)*[x-xs, y-ys, z-zs];
rangeRatePartialV = [x-xs, y-ys, z-zs]/range;

% station location partials
% rangePartialRs = -[x-xs, y-ys, z-zs]/range;
rangePartialRs = [x-xs, y-ys, z-zs]/range*ECEF2ECI;
% rangeRatePartialRs = -[xdot-xsdot, ydot-ysdot, zdot-zsdot]/range + rangeRate/(range^2)*[x-xs, y-ys, z-zs];
rangeRatePartialRs = (-[xdot-xsdot, ydot-ysdot, zdot-zsdot]/range + rangeRate/(range^2)*[x-xs, y-ys, z-zs])*ECEF2ECI;


H_range = [rangePartialR, rangePartialV, rangePartialRs];
H_rangeRate = [rangeRatePartialR, rangeRatePartialV, rangeRatePartialRs];
end

