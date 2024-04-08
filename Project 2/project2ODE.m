function [Sdot] = project2ODE(t,S,const)
%ORBITODE Summary of this function goes here
%   S is the 6x1 state vector [position; velocity]

muE = const.muE;
muS = const.muS;
ae = const.ae;

n = 7;

x = S(1);
y = S(2);
z = S(3);
rE_N = S(1:3);
rE = norm(rE_N);
v_N = S(4:6);
Cr = S(7);

Phi = reshape(S(n+1:end),[n,n]);

% Dynamics
aE_N = -muE*rE_N/(rE^3);

JD = const.t0 + t/86400; % JD (t must be in seconds)
[r_S2E_N, ~, ~] = Ephem(JD,3,'EME2000');
r_S2E = norm(r_S2E_N);
rS_N = r_S2E_N + rE_N; % position vector from sun to spacecraft
rS = norm(rS_N);
% aS_N = -muS*rS_N/(rS^3);
aS_N = muS*(-rS_N/rS^3 + r_S2E_N/r_S2E^3);

p_sr = const.p_srAU*(const.AU^2)/(rS^2);
aSRP_N = Cr*p_sr*const.a2m/rS*rS_N;

a_N = aE_N + aS_N + aSRP_N;

% Jacobians
amuPartialR = [
    -((sqrt(rE^2) * (rE^2 - 3*x^2) * muE) / rE^6), (3*x*y*muE) / (rE^2)^(5/2), (3*x*z*muE) / (rE^2)^(5/2);
    (3*x*y*muE) / (rE^2)^(5/2), -((sqrt(rE^2) * (rE^2 - 3*y^2) * muE) / rE^6), (3*y*z*muE) / (rE^2)^(5/2);
    (3*x*z*muE) / (rE^2)^(5/2), (3*y*z*muE) / (rE^2)^(5/2), -((sqrt(rE^2) * (rE^2 - 3*z^2) * muE) / rE^6)
];

% A = [
%     zeros(3), eye(3);
%     amuPartialR, zeros(3);
% ];
A = zeros(n);


PhiDot = A*Phi;

Sdot = [v_N; a_N; 0; reshape(PhiDot,[],1)];
end

