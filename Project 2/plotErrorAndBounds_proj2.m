function [] = plotErrorAndBounds_proj2(t, StateErrors, P, title)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_1 = squeeze(P(1,1,:)).^(1/2);
sig_2 = squeeze(P(2,2,:)).^(1/2);
sig_3 = squeeze(P(3,3,:)).^(1/2);
sig_4 = squeeze(P(4,4,:)).^(1/2);
sig_5 = squeeze(P(5,5,:)).^(1/2);
sig_6 = squeeze(P(6,6,:)).^(1/2);
sig_7 = squeeze(P(7,7,:)).^(1/2);

maxError = max(abs(StateErrors),[], 2);

t = t./60^2;

figure
subplot(4, 2, 1);
sgtitle(title, 'Interpreter', 'latex');
% plot(t, StateErrors(1, :))
scatter(t, StateErrors(1, :), '.')
hold on
plot(t, 3*sig_1, 'r--')
plot(t, -3*sig_1, 'r--')
xlabel('time [hours]');
ylabel("$\delta x$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(1))
% legend("STM", "integrated", 'Location', 'northwest')
% subplot
subplot(4, 2, 3);
% plot(t, StateErrors(2, :));
scatter(t, StateErrors(2, :), '.');
hold on
plot(t, 3*sig_2, 'r--')
plot(t, -3*sig_2, 'r--')
xlabel('time [hours]');
ylabel("$\delta y$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(2))
% subplot
subplot(4, 2, 5);
% plot(t, StateErrors(3, :));
scatter(t, StateErrors(3, :), '.');
hold on
plot(t, 3*sig_3, 'r--')
plot(t, -3*sig_3, 'r--')
xlabel('time [hours]');
ylabel("$\delta z$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(3))
% subplot
subplot(4, 2, 2);
% plot(t, StateErrors(4, :));
scatter(t, StateErrors(4, :), '.');
hold on
plot(t, 3*sig_4, 'r--')
plot(t, -3*sig_4, 'r--')
xlabel('time [hours]');
ylabel("$\delta \dot{x}$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(4))
% subplot
subplot(4, 2, 4);
% plot(t, StateErrors(5, :));
scatter(t, StateErrors(5, :), '.');
hold on
plot(t, 3*sig_5, 'r--')
plot(t, -3*sig_5, 'r--')
xlabel('time [hours]');
ylabel("$\delta \dot{y}$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(5))
% subplot
subplot(4, 2, 6);
% plot(t, StateErrors(6, :));
scatter(t, StateErrors(6, :), '.');
hold on
plot(t, 3*sig_6, 'r--')
plot(t, -3*sig_6, 'r--')
xlabel('time [hours]');
ylabel("$\delta \dot{z}$", 'Interpreter', 'latex');
ylim([-1 1]*maxError(6))
subplot(4, 2, 7);
% plot(t, StateErrors(6, :));
scatter(t, StateErrors(7, :), '.');
hold on
plot(t, 3*sig_7, 'r--')
plot(t, -3*sig_7, 'r--')
xlabel('time [hours]');
ylabel("$\delta C_R$", 'Interpreter', 'latex');
% ylim([-1 1]*maxError(7))
end

