function [] = plotResiduals_HW3(t, y, alpha, noise_sd, titles)
%PLOTRESIDUALS Summary of this function goes here
%   Detailed explanation goes here

t = t./60^2;

figure
subplot(2,1,1)
sgtitle(titles(1), 'Interpreter', 'latex')
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
xlabel("Time [hours]")
% legend("$\rho$", "$\dot{\rho}$", 'Interpreter', 'latex')

figure
subplot(2,1,1)
sgtitle(titles(2), 'Interpreter', 'latex')
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
xlabel("Time [hours]")
end

