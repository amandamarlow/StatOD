function [omegaRcN_N] = omegaRcN_N(t)
addpath("C:\Users\marlo\MATLAB Drive\6010\RigidBodyKinematics-Matlab\Matlab")
    % numerical method to get RcNdot
    deltat = 1*10^(-8);
    t0 = t-deltat/2;
    t1 = t+deltat/2;
    RcNdot = (RcN(t1)-RcN(t0))/deltat;
    %find omegaRcN_N
    omegatilde = -RcN(t)'*RcNdot;
    omegaRcN_N = [-omegatilde(2,3); omegatilde(1,3); -omegatilde(1,2)];
end
