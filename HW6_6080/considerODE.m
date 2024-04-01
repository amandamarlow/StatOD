function [Sdot] = considerODE(t,S,n,q,constants)
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

% mu = S(7);
% J2 = S(8);
% J3 = S(9);
% C = [mu; J2; J3];

Phi = reshape(S(n+1:n+n^2),[n,n]);
theta = reshape(S(end+1-n*q:end),[n,q]);

% Dynamics
a_mu_N = -mu*r_N/(r^3);
a_J2_N = -3*mu*(ae^2)*J2/2/(r^7) * [x*(r^2 - 5*z^2); y*(r^2 - 5*z^2); z*(r^2 + 2*(x^2 + y^2) - 3*z^2)];
a_N = a_mu_N + a_J2_N;

% Jacobians
amuPartialR = [
    -((sqrt(r^2) * (r^2 - 3*x^2) * mu) / r^6), (3*x*y*mu) / (r^2)^(5/2), (3*x*z*mu) / (r^2)^(5/2);
    (3*x*y*mu) / (r^2)^(5/2), -((sqrt(r^2) * (r^2 - 3*y^2) * mu) / r^6), (3*y*z*mu) / (r^2)^(5/2);
    (3*x*z*mu) / (r^2)^(5/2), (3*y*z*mu) / (r^2)^(5/2), -((sqrt(r^2) * (r^2 - 3*z^2) * mu) / r^6)
];
aJ2PartialR = [
    -((3 * J2 * ae^2 * (-7*x^2*(x^2+y^2-4*z^2) + r^2*(3*x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2))), (3 * J2 * ae^2 * x*y * (-2*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * x*z * (8*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2));
    (3 * J2 * ae^2 * x*y * (-2*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2)), -((3 * J2 * ae^2 * (-7*y^2*(x^2+y^2-4*z^2) + r^2*(x^2+3*y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2))), (3 * J2 * ae^2 * y*z * (8*r^2 + 7*(x^2+y^2-4*z^2)) * mu) / (2 * r^8 * sqrt(r^2));
    (3 * J2 * ae^2 * x*z * (-6*r^2 + 21*(x^2+y^2) - 14*z^2) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * y*z * (-6*r^2 + 21*(x^2+y^2) - 14*z^2) * mu) / (2 * r^8 * sqrt(r^2)), (3 * J2 * ae^2 * (21*(x^2+y^2)*z^2 - 14*z^4 - 3*r^2*(x^2+y^2-2*z^2)) * mu) / (2 * r^8 * sqrt(r^2))
];


% amuPartialC = [
%     -(x/r^3), 0, 0;
%     -(y/r^3), 0, 0;
%     -(z/r^3), 0, 0
% ];
% aJ2PartialC = [
%     -((3 * J2 * ae^2 * x * (x^2+y^2-4*z^2)) / (2 * (x^2+y^2+z^2)^(7/2))), -((3 * ae^2 * x * (x^2+y^2-4*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2))), 0;
%     -((3 * J2 * ae^2 * y * (x^2+y^2-4*z^2)) / (2 * (x^2+y^2+z^2)^(7/2))), -((3 * ae^2 * y * (x^2+y^2-4*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2))), 0;
%     (3 * J2 * ae^2 * z * (-3*(x^2+y^2)+2*z^2)) / (2 * (x^2+y^2+z^2)^(7/2)), (3 * ae^2 * z * (-3*(x^2+y^2)+2*z^2) * mu) / (2 * (x^2+y^2+z^2)^(7/2)), 0
% ];


aJ3PartialJ3 = 1/2*mu/r^2*(ae/r)^3 * [5*(7*(z/r)^3 - 3*(z/r))*x/r; 5*(7*(z/r)^3 - 3*(z/r))*y/r; 3*(1 - 10*(z/r)^2 + 35/3*(z/r)^4)];

% A = [
%     zeros(3), eye(3), zeros(3);
%     amuPartialR+aJ2PartialR, zeros(3), amuPartialC+aJ2PartialC;
%     zeros(3), zeros(3), zeros(3);
% ];
A = [
    zeros(3), eye(3);
    amuPartialR+aJ2PartialR, zeros(3);
];

% Ac = amuPartialC + aJ2PartialC;
B = [zeros(3, 1); aJ3PartialJ3];

PhiDot = A*Phi;
thetaDot = A*theta + B;

Sdot = [v_N; a_N; reshape(PhiDot,[],1); reshape(thetaDot,[],1)];
end

