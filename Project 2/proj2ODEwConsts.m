function [Sdot] = proj2ODEwConsts(t,S,const)
%ORBITODE Summary of this function goes here
%   S is the 6x1 state vector [position; velocity]

% muE = const.muE;
% muS = const.muS;
p_srAU = const.p_srAU;
a2m = const.a2m;
% ae = const.ae;

% q = 14; % # of constants
% q = 13; % # of constants
q = 12; % # of constants
% q = 11; % # of constants
n = 6+q;

x = S(1);
y = S(2);
z = S(3);
rE_N = S(1:3);
rE = norm(rE_N);
v_N = S(4:6);
Cr = S(7);
muE = S(8);
muS = S(9);
% p_srAU = S(10);
% a2m = S(11);



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

p_sr = p_srAU*(const.AU^2)/(r_S2sc^2);
aSRP_N = Cr*p_sr*a2m/r_S2sc*r_S2sc_N;

aE_N = -muE*rE_N/(rE^3);

a_N = aE_N + aS_N + aSRP_N;

% Jacobians
amu_partial_R = -muE/(rE^3)*eye(3) + 3*muE/(rE^5)*(rE_N*rE_N');
aS_partial_R = muS * (3/(r_S2sc^5)*(r_sc2S_N*r_sc2S_N')*eye(3) - 1/r_S2sc^3*eye(3));
aSRP_partial_R = Cr*p_srAU*(const.AU^2)*a2m * (1/(r_S2sc^3)*eye(3) - 3/(r_S2sc^5)*(r_S2sc_N*r_S2sc_N')*eye(3));

aSRP_partial_Cr = p_sr*a2m/r_S2sc*r_S2sc_N;
aE_partial_muE = -rE_N/(rE^3);
aS_partial_muS = (-r_S2sc_N/r_S2sc^3 + r_S2E_N/r_S2E^3);
% aSRP_partial_p = Cr*(const.AU^2)/(r_S2sc^3)*a2m*r_S2sc_N;
% aSRP_partial_A2m = Cr*p_sr/r_S2sc*r_S2sc_N;

a_partial_R = amu_partial_R + aS_partial_R + aSRP_partial_R;

% a_partial_c = [aSRP_partial_Cr, aE_partial_muE, aS_partial_muS, aSRP_partial_p, aSRP_partial_A2m, zeros(3,9)];
% a_partial_c = [aSRP_partial_Cr, aE_partial_muE, aS_partial_muS, aSRP_partial_p, zeros(3,9)];
a_partial_c = [aSRP_partial_Cr, aE_partial_muE, aS_partial_muS, zeros(3,9)];
% a_partial_c = [aSRP_partial_Cr, aE_partial_muE, zeros(3,9)];
A = [
    zeros(3), eye(3), zeros(3,q);
    a_partial_R, zeros(3), a_partial_c;
    zeros(q,n)
];


PhiDot = A*Phi;

Sdot = [v_N; a_N; zeros(q,1); reshape(PhiDot,[],1)];
end

