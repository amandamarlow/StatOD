function [r_E] = latlon2ECEF(r, phi, lambda)
%LATLON2ECEF Summary of this function goes here
%   Detailed explanation goes here

x = r*cos(phi)*cos(lambda);
y = r*cos(phi)*sin(lambda);
z = r*sin(phi);

r_E = [x; y; z];
end
