function [] = plotErrorAndBounds(StateErrors,P)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_rng = squeeze(P(1,1,:)).^(1/2);
sig_rngRt = squeeze(P(2,2,:)).^(1/2);

figure
subplot(3, 2, 1);
sgtitle('Deviation vectors vs time');
plot(t2, StateErrors(1, :))
hold on
yline(3*sig_rng(1), 'r--')
yline(-3*sig_rng(1), 'r--')
xlabel('time [s]');
ylabel("$\delta x$", 'Interpreter', 'latex');
legend("STM", "integrated", 'Location', 'northwest')
% subplot
subplot(3, 2, 3);
plot(t2, StateErrors(2, :));
hold on
yline(3*sig_rng(1), 'r--')
yline(-3*sig_rng(1), 'r--')
xlabel('time [s]');
ylabel("$\delta y$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 5);
plot(t2, StateErrors(3, :));
hold on
yline(3*sig_rng(1), 'r--')
yline(-3*sig_rng(1), 'r--')
xlabel('time [s]');
ylabel("$\delta z$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 2);
plot(t2, StateErrors(4, :));
hold on
yline(3*sig_rngRt(1), 'r--')
yline(-3*sig_rngRt(1), 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{x}$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 4);
plot(t2, StateErrors(5, :));
hold on
yline(3*sig_rngRt(1), 'r--')
yline(-3*sig_rngRt(1), 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{y}$", 'Interpreter', 'latex');
% subplot
subplot(3, 2, 6);
plot(t2, StateErrors(6, :));
hold on
yline(3*sig_rngRt(1), 'r--')
yline(-3*sig_rngRt(1), 'r--')
xlabel('time [s]');
ylabel("$\delta \dot{z}$", 'Interpreter', 'latex');
end

