function [] = plotJ3aError(t, StateErrors, P, title)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_1 = squeeze(P(1,1,:)).^(1/2);
sig_2 = squeeze(P(2,2,:)).^(1/2);
sig_3 = squeeze(P(3,3,:)).^(1/2);

t = t./60^2;

figure
subplot(3, 1, 1);
sgtitle(title, 'Interpreter', 'latex');
% plot(t, StateErrors(1, :))
scatter(t, StateErrors(1, :), '.')
hold on
plot(t, 3*sig_1, 'r--')
plot(t, -3*sig_1, 'r--')
xlabel('time [hours]');
ylabel("$\delta x$", 'Interpreter', 'latex');
% ylim([-1 1]*maxError(1))
% legend("STM", "integrated", 'Location', 'northwest')
% subplot
subplot(3, 1, 2);
% plot(t, StateErrors(2, :));
scatter(t, StateErrors(2, :), '.');
hold on
plot(t, 3*sig_2, 'r--')
plot(t, -3*sig_2, 'r--')
xlabel('time [hours]');
ylabel("$\delta y$", 'Interpreter', 'latex');
% ylim([-1 1]*maxError(2))
% subplot
subplot(3, 1, 3);
% plot(t, StateErrors(3, :));
scatter(t, StateErrors(3, :), '.');
hold on
plot(t, 3*sig_3, 'r--')
plot(t, -3*sig_3, 'r--')
xlabel('time [hours]');
ylabel("$\delta z$", 'Interpreter', 'latex');
end

