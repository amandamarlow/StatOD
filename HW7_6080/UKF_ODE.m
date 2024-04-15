function [Sdot] = UKF_ODE(t,S,n,constants)
%ORBITODE Summary of this function goes here
%   S is the 6x1 state vector [position; velocity]

mu = constants.mu;
ae = constants.ae;
J2 = constants.J2;

x = S(1);
y = S(2);
z = S(3);
r_N = S(1:3);
r = norm(r_N);
v_N = S(4:6);

% Dynamics
a_mu_N = -mu*r_N/(r^3);
a_J2_N = -3*mu*(ae^2)*J2/2/(r^7) * [x*(r^2 - 5*z^2); y*(r^2 - 5*z^2); z*(r^2 + 2*(x^2 + y^2) - 3*z^2)];
a_N = a_mu_N + a_J2_N;


Sdot = [v_N; a_N];
end

