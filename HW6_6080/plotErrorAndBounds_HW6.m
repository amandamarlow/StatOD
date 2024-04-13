function [] = plotErrorAndBounds_HW6(t, StateErrors, P, Pxx, title)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_1 = squeeze(P(1,1,:)).^(1/2);
sig_2 = squeeze(P(2,2,:)).^(1/2);
sig_3 = squeeze(P(3,3,:)).^(1/2);
sig_4 = squeeze(P(4,4,:)).^(1/2);
sig_5 = squeeze(P(5,5,:)).^(1/2);
sig_6 = squeeze(P(6,6,:)).^(1/2);

sig_1xx = squeeze(Pxx(1,1,:)).^(1/2);
sig_2xx = squeeze(Pxx(2,2,:)).^(1/2);
sig_3xx = squeeze(Pxx(3,3,:)).^(1/2);
sig_4xx = squeeze(Pxx(4,4,:)).^(1/2);
sig_5xx = squeeze(Pxx(5,5,:)).^(1/2);
sig_6xx = squeeze(Pxx(6,6,:)).^(1/2);

t = t./60^2;

maxError = max(abs(StateErrors),[], 2);

% %colororder("reef");

figure
subplot(3, 2, 1);
sgtitle(title, 'Interpreter', 'latex');
% plot(t, StateErrors(1, :))
scatter(t, StateErrors(1, :), '.')
hold on
Pplus2 = plot(t, 2*sig_1, '--');
Pminus2 = plot(t, -2*sig_1, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 2*sig_1xx, '--');
Pminus2 = plot(t, -2*sig_1xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta x$", 'Interpreter', 'latex');
%colororder("reef");
ylim([-1 1]*maxError(1))
% subplot
subplot(3, 2, 3);
% plot(t, StateErrors(2, :));
scatter(t, StateErrors(2, :), '.');
hold on
Pplus2 = plot(t, 2*sig_2, '--');
Pminus2 = plot(t, -2*sig_2, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 2*sig_2xx, '--');
Pminus2 = plot(t, -2*sig_2xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta y$", 'Interpreter', 'latex');
%colororder("reef");
ylim([-1 1]*maxError(2))
% subplot
subplot(3, 2, 5);
% plot(t, StateErrors(3, :));
scatter(t, StateErrors(3, :), '.');
hold on
Pplus2 = plot(t, 2*sig_3, '--');
Pminus2 = plot(t, -2*sig_3, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 2*sig_3xx, '--');
Pminus2 = plot(t, -2*sig_3xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta z$", 'Interpreter', 'latex');
%colororder("reef");
ylim([-1 1]*maxError(3))
% subplot
subplot(3, 2, 2);
% plot(t, StateErrors(4, :));
scatter(t, StateErrors(4, :), '.');
hold on
Pplus2 = plot(t, 2*sig_4, '--');
Pminus2 = plot(t, -2*sig_4, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 2*sig_4xx, '--');
Pminus2 = plot(t, -2*sig_4xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{x}$", 'Interpreter', 'latex');
%colororder("reef");
legend("$\delta x$", "2$\sigma$ ($P_x$)", '', "2$\sigma$ ($P_{xx}$)", '','Location', 'northeast', "Interpreter", "latex")
ylim([-1 1]*maxError(4))
% subplot
subplot(3, 2, 4);
% plot(t, StateErrors(5, :));
scatter(t, StateErrors(5, :), '.');
hold on
Pplus1 = plot(t, 2*sig_5, '--');
Pminus1 = plot(t, -2*sig_5, '--');
Pminus1.SeriesIndex = Pplus1.SeriesIndex;
Pplus2 = plot(t, 2*sig_5xx, '--');
Pminus2 = plot(t, -2*sig_5xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{y}$", 'Interpreter', 'latex');
%colororder("reef");
ylim([-1 1]*maxError(5))
% subplot
subplot(3, 2, 6);
% plot(t, StateErrors(6, :));
scatter(t, StateErrors(6, :), '.');
hold on
Pplus2 = plot(t, 2*sig_6, '--');
Pminus2 = plot(t, -2*sig_6, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 2*sig_6xx, '--');
Pminus2 = plot(t, -2*sig_6xx, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{z}$", 'Interpreter', 'latex');
%colororder("reef");
ylim([-1 1]*maxError(6))
end
