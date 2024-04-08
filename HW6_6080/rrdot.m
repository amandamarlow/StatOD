function [r_N, rdot_N] = rrdot(X)
% 
mu = 42828.3; % [km^3/s^2] mars gravitational constant

r_mag = X(1);
OMEGA = X(2); 
i = X(3);
theta = X(4);

r_O = r_mag*[1;0;0]; % In hill/orbit frame
thetadot = sqrt(mu/r_mag^3);
omegaON = thetadot*[0;0;1]; % Angular velocity of H relative to N written in O coordinates

ON = Euler3132C([OMEGA; i; theta]);
NO = ON';

r_N = NO*r_O; % position vector in inertial frame
% use transport theorem for velocity. Note: rdot_H is assumed 0 for this circular orbit
rdot_N = NO*cross(omegaON,r_O); % velocity vector in inertial frame
end
