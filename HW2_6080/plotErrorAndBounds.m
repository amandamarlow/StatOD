function [] = plotErrorAndBounds(t, StateErrors, P, title)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_1 = squeeze(P(1,1,:)).^(1/2);
sig_2 = squeeze(P(2,2,:)).^(1/2);
sig_3 = squeeze(P(3,3,:)).^(1/2);
sig_4 = squeeze(P(4,4,:)).^(1/2);
sig_5 = squeeze(P(5,5,:)).^(1/2);
sig_6 = squeeze(P(6,6,:)).^(1/2);

figure
subplot(3, 2, 1);
sgtitle(title, 'Interpreter', 'latex');
% plot(t, StateErrors(1, :))
scatter(t, StateErrors(1, :), '.')
hold on
plot(t, 3*sig_1, 'r--')
plot(t, -3*sig_1, 'r--')
xlabel('time [s]');
ylabel("$\delta x$", 'Interpreter', 'latex');
% legend("STM", "integrated", 'Location', 'northwest')
% subplot
subplot(3, 2, 3);
% plot(t, StateErrors(2, :));
scatter(t, StateErrors(2, :), '.');
hold on
plot(t, 3*sig_2, 'r--')
plot(t, -3*sig_2, 'r--')
xlabel('time [s]');
ylabel("$\delta y$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 5);
% plot(t, StateErrors(3, :));
scatter(t, StateErrors(3, :), '.');
hold on
plot(t, 3*sig_3, 'r--')
plot(t, -3*sig_3, 'r--')
xlabel('time [s]');
ylabel("$\delta z$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 2);
% plot(t, StateErrors(4, :));
scatter(t, StateErrors(4, :), '.');
hold on
plot(t, 3*sig_4, 'r--')
plot(t, -3*sig_4, 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{x}$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 4);
% plot(t, StateErrors(5, :));
scatter(t, StateErrors(5, :), '.');
hold on
plot(t, 3*sig_5, 'r--')
plot(t, -3*sig_5, 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{y}$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 6);
% plot(t, StateErrors(6, :));
scatter(t, StateErrors(6, :), '.');
hold on
plot(t, 3*sig_6, 'r--')
plot(t, -3*sig_6, 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{z}$", 'Interpreter', 'latex');
end

