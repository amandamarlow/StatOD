function [Sdot] = hw8ODE(t,S)

X = S(1:2);
Phi = reshape(S(3:end), 2,2);
A = [0, 1; 0, 0];

Xdot = A*X;
PhiDot = A*Phi;

Sdot = [Xdot, reshape(PhiDot, [], 1)];
end

