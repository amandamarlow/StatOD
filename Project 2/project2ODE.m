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

JD = const.t0 + t/86400; % JD (t must be in seconds)
[r_S2E_N, ~, ~] = Ephem(JD,3,'EME2000');
r_S2E = norm(r_S2E_N);
r_S2sc_N = r_S2E_N + rE_N; % position vector from sun to spacecraft
r_sc2S_N = -r_S2sc_N;
r_S2sc = norm(r_S2sc_N);
% aS_N = -muS*rS_N/(rS^3);
aS_N = muS*(-r_S2sc_N/r_S2sc^3 + r_S2E_N/r_S2E^3);

p_sr = const.p_srAU*(const.AU^2)/(r_S2sc^2);
aSRP_N = Cr*p_sr*const.a2m/r_S2sc*r_S2sc_N;

aE_N = -muE*rE_N/(rE^3);

a_N = aE_N + aS_N + aSRP_N;

% Jacobians
amu_partial_R = [
    -((sqrt(rE^2) * (rE^2 - 3*x^2) * muE) / rE^6), (3*x*y*muE) / (rE^2)^(5/2), (3*x*z*muE) / (rE^2)^(5/2);
    (3*x*y*muE) / (rE^2)^(5/2), -((sqrt(rE^2) * (rE^2 - 3*y^2) * muE) / rE^6), (3*y*z*muE) / (rE^2)^(5/2);
    (3*x*z*muE) / (rE^2)^(5/2), (3*y*z*muE) / (rE^2)^(5/2), -((sqrt(rE^2) * (rE^2 - 3*z^2) * muE) / rE^6)
];
aS_partial_R = muS * (3/(r_S2sc^5)*(r_sc2S_N*r_sc2S_N')*eye(3) - 1/r_S2sc^3*eye(3));
aSRP_partial_R = Cr*p_sr*const.a2m * (1/(r_S2sc^3)*eye(3) - 3/(r_S2sc^5)*(r_S2sc_N*r_S2sc_N')*eye(3));

aSRP_partial_Cr = p_sr*const.a2m/r_S2sc*r_S2sc_N;

a_partial_R = amu_partial_R + aS_partial_R + aSRP_partial_R;
a_partial_c = aSRP_partial_Cr;
A = [
    zeros(3), eye(3), zeros(3,1);
    a_partial_R, zeros(3), a_partial_c;
    zeros(1,n)
];


PhiDot = A*Phi;

Sdot = [v_N; a_N; 0; reshape(PhiDot,[],1)];
end

