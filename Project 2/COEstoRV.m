function [R_N,V_N] = COEstoRV(a,e,i,OMEGA,omega,nu,mu_s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = a*(1-e^2);
R = p/(1+e*cos(nu));
R_hill = [R; 0; 0];
h = sqrt(mu_s*p);
V_hill = [mu_s/h*e*sin(nu); mu_s/h*(1+e*cos(nu)); 0];
PN = M3(omega+nu)*M1(i)*M3(OMEGA);
R_N = PN'*R_hill; 
V_N = PN'*V_hill;
end

