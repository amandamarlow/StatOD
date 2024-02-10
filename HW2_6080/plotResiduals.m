function [] = plotResiduals(t, y, alpha, noise_sd, titles)
%PLOTRESIDUALS Summary of this function goes here
%   Detailed explanation goes here

figure
subplot(2,1,1)
sgtitle(titles(1))
scatter(t, y(1,:), '.')
hold on
yline(3*noise_sd(1), 'r--')
yline(-3*noise_sd(1), 'r--')
ylabel("$\rho$ residual [km]", 'Interpreter', 'latex')
subplot(2,1,2)
scatter(t, y(2,:), '.')
hold on
yline(3*noise_sd(2), 'r--')
yline(-3*noise_sd(2), 'r--')
ylabel("$\dot{\rho}$ residual [km/s]", 'Interpreter', 'latex')
xlabel("Time [s]")
% legend("$\rho$", "$\dot{\rho}$", 'Interpreter', 'latex')

figure
subplot(2,1,1)
sgtitle(titles(2))
scatter(t, alpha(1,:), '.')
hold on
yline(3*noise_sd(1), 'r--')
yline(-3*noise_sd(1), 'r--')
ylabel("$\rho$ residual [km]", 'Interpreter', 'latex')
subplot(2,1,2)
scatter(t, alpha(2,:), '.')
hold on
yline(3*noise_sd(2), 'r--')
yline(-3*noise_sd(2), 'r--')
ylabel("$\dot{\rho}$ residual [km/s]", 'Interpreter', 'latex')
xlabel("Time [s]")
end

