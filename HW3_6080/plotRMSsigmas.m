function [] = plotRMSsigmas(sigmas, RMSpos, RMSvel, RMSresiduals, title)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
subplot(2,1,1)
sgtitle(title)
% semilogx(sigmas, RMSpos, 'o-')
loglog(sigmas, RMSpos, 'o-')
hold on
% semilogx(sigmas, RMSresiduals(1,:), 'o-')
loglog(sigmas, RMSresiduals(1,:), 'o-')
legend(["Position", "Range"])
xlabel("$\sigma$ $[km]$", 'Interpreter','latex')
ylabel("RMS [km]")
subplot(2,1,2)
% semilogx(sigmas, RMSvel, 'o-')
loglog(sigmas, RMSvel, 'o-')
hold on
% semilogx(sigmas, RMSresiduals(2,:), 'o-')
loglog(sigmas, RMSresiduals(2,:), 'o-')
legend("Velocity", "Range Rate")
xlabel("$\sigma$ $[km/s^2]$", 'Interpreter','latex')
ylabel("RMS [km/s]")
end

