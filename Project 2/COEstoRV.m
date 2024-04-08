function [R,V] = COEstoRV(a,e,i,OMEGA,omega,nu,mu_s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rp = a*(1-e);
R_perifocal = [rp; 0; 0];
vp = sqrt((1+e)/(1-e)*mu_s/a);
v0_perifocal = [0; vp; 0];
PN = M3(omega+nu)*M1(i)*M3(OMEGA);
R = PN'*R_perifocal; 
V = PN'*v0_perifocal;
end

