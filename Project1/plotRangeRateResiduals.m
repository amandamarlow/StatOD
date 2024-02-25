function [] = plotRangeRateResiduals(t, y, alpha, noise_sd, titles)
%PLOTRESIDUALS Summary of this function goes here
%   Detailed explanation goes here

t = t./60^2;

% figure
% scatter(t, y, '.')
% hold on
% yline(3*noise_sd(2), 'r--')
% yline(-3*noise_sd(2), 'r--')
% ylabel("$\dot{\rho}$ residual [km/s]", 'Interpreter', 'latex')
% xlabel("Time [hours]")
% title(titles(1), 'Interpreter', 'latex')

figure
scatter(t, alpha, '.')
hold on
yline(3*noise_sd(2), 'r--')
yline(-3*noise_sd(2), 'r--')
ylabel("$\dot{\rho}$ residual [km/s]", 'Interpreter', 'latex')
xlabel("Time [hours]")
title(titles(2), 'Interpreter', 'latex')
end

