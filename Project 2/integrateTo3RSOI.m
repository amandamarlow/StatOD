function [t, Xvec, STM_vec] = integrateTo3RSOI(tspan, X_t1, STM_t1, constants)
% function [t, Xvec] = integrateTo3RSOI(tspan, X_t1, STM_t1, constants)
%INTEGRATETRAJECTORYSTM Summary of this function goes here
%   Detailed explanation goes here

n = length(X_t1);
% tspan = [t1 t2];
S0 = [X_t1; reshape(STM_t1,[],1)];

% simulate of reference trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @(t,S) RSOIeventsFcn(t,S,constants));
[t,S] = ode45(@(t,S) project2ODE(t,S,constants), tspan, S0, options);
% [t,S] = ode45(@(t,S) rvOnlyODE(t,S,constants), tspan, S0, options);

% STM
flatSTM = S(:,n+1:end);
STM_vec = reshape(flatSTM',n,n,[]);

% return
Xvec = S(:,1:n)';
t = t';
end

function [position,isterminal,direction] = RSOIeventsFcn(t,y,constants)
  RSOI = constants.RSOI;
  position = norm(y(1:3)) - 3*RSOI; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

% function [Sdot]= rvOnlyODE(t,S,const)
% %ORBITODE Summary of this function goes here
% %   S is the 6x1 state vector [position; velocity]
% 
% % muE = const.muE;
% muS = const.muS;
% ae = const.ae;
% 
% n = 7;
% 
% x = S(1);
% y = S(2);
% z = S(3);
% rE_N = S(1:3);
% rE = norm(rE_N);
% v_N = S(4:6);
% Cr = S(7);
% 
% % Dynamics
% 
% JD = const.t0 + t/86400; % JD (t must be in seconds)
% [r_S2E_N, ~, muE] = Ephem(JD,3,'EME2000');
% r_S2E = norm(r_S2E_N);
% r_S2sc_N = r_S2E_N + rE_N; % position vector from sun to spacecraft
% r_sc2S_N = -r_S2sc_N;
% r_S2sc = norm(r_S2sc_N);
% % aS_N = -muS*rS_N/(rS^3);
% aS_N = muS*(-r_S2sc_N/r_S2sc^3 + r_S2E_N/r_S2E^3);
% 
% p_sr = const.p_srAU*(const.AU^2)/(r_S2sc^2);
% aSRP_N = Cr*p_sr*const.a2m/r_S2sc*r_S2sc_N;
% aE_N = -muE*rE_N/(rE^3);
% 
% a_N = aE_N + aS_N + aSRP_N;
% 
% Sdot = [v_N; a_N];
% end
