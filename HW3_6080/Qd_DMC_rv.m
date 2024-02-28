function [Qd] = Qd_DMC_rv(tau, dt, Qc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pp = tau^5/2 * ((1 - exp(-2*dt/tau)) + 2*dt/tau*(1 - 2*exp(-dt/tau)) - 2*(dt/tau)^2 + 2/3*(dt/tau)^3);
pv = tau^4/2 * ((exp(-2*dt/tau) - 1) - 2*(exp(-dt/tau) - 1) + 2*dt/tau*(exp(-dt/tau) - 1) + (dt/tau)^2);
pa = tau^3/2 * ((1 - exp(-2*dt/tau)) - 2*dt/tau*exp(-dt/tau));
vv = tau^3/2 * ((1-exp(-2*dt*tau)) - 4*(1 - exp(-dt/tau)) + 2*dt/tau);
va = tau^2/2 * (1 - exp(-dt/tau))^2;
aa = tau/2 * (1 - exp(-2*dt/tau));

Qd = [pp*Qc, pv*Qc, pa*Qc;
    pv*Qc, vv*Qc, va*Qc;
    pa*Qc, va*Qc, aa*Qc];

end