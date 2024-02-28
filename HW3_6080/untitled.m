function [outputArg1,outputArg2] = plotRMSsigmas()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
subplot(2,1,1)
sgtitle("CKF-SNC RMS Values")
semilogx(sigmas, sigsRMSposition_CKF_SNC, 'o-')
hold on
semilogx(sigmas, sigsRMSresiduals_CKF_SNC(1,:), 'o-')
legend(["Position", "Range"])
xlabel("$\sigma$ $[km]$", 'Interpreter','latex')
ylabel("RMS [km]")
subplot(2,1,2)
semilogx(sigmas, sigsRMSvelocity_CKF_SNC, 'o-')
hold on
semilogx(sigmas, sigsRMSresiduals_CKF_SNC(2,:), 'o-')
legend("Velocity", "Range Rate")
xlabel("$\sigma$ $[km/s^2]$", 'Interpreter','latex')
ylabel("RMS [km/s]")
end

