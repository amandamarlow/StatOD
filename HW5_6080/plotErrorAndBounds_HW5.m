function [] = plotErrorAndBounds_HW5(t, StateErrors_CKF, P_CKF, StateErrors_batch, P_batch, title)
%PLOTERRORANDBOUNDS Summary of this function goes here
%   Detailed explanation goes here

sig_1_CKF = squeeze(P_CKF(1,1,:)).^(1/2);
sig_2_CKF = squeeze(P_CKF(2,2,:)).^(1/2);
sig_3_CKF = squeeze(P_CKF(3,3,:)).^(1/2);
sig_4_CKF = squeeze(P_CKF(4,4,:)).^(1/2);
sig_5_CKF = squeeze(P_CKF(5,5,:)).^(1/2);
sig_6_CKF = squeeze(P_CKF(6,6,:)).^(1/2);

sig_1_batch = squeeze(P_batch(1,1,:)).^(1/2);
sig_2_batch = squeeze(P_batch(2,2,:)).^(1/2);
sig_3_batch = squeeze(P_batch(3,3,:)).^(1/2);
sig_4_batch = squeeze(P_batch(4,4,:)).^(1/2);
sig_5_batch = squeeze(P_batch(5,5,:)).^(1/2);
sig_6_batch = squeeze(P_batch(6,6,:)).^(1/2);

t = t./60^2;

% colororder("reef");

figure
subplot(3, 2, 1);
sgtitle(title, 'Interpreter', 'latex');
% plot(t, StateErrors(1, :))
scatter(t, StateErrors_CKF(1, :), 'x')
hold on
scatter(t, StateErrors_batch(1, :), '.')
Pplus2 = plot(t, 3*sig_1_CKF, '--');
Pminus2 = plot(t, -3*sig_1_CKF, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 3*sig_1_batch, '--');
Pminus2 = plot(t, -3*sig_1_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta x$", 'Interpreter', 'latex');
colororder("reef");
% ylim([-1 1]*maxError(1))
% subplot
subplot(3, 2, 3);
% plot(t, StateErrors(2, :));
scatter(t, StateErrors_CKF(2, :), 'x');
hold on
scatter(t, StateErrors_batch(2, :), '.');
Pplus2 = plot(t, 3*sig_2_CKF, '--');
Pminus2 = plot(t, -3*sig_2_CKF, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 3*sig_2_batch, '--');
Pminus2 = plot(t, -3*sig_2_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta y$", 'Interpreter', 'latex');
colororder("reef");
% ylim([-1 1]*maxError(2))
% subplot
subplot(3, 2, 5);
% plot(t, StateErrors(3, :));
scatter(t, StateErrors_CKF(3, :), 'x');
hold on
scatter(t, StateErrors_batch(3, :), '.');
Pplus2 = plot(t, 3*sig_3_CKF, '--');
Pminus2 = plot(t, -3*sig_3_CKF, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 3*sig_3_batch, '--');
Pminus2 = plot(t, -3*sig_3_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta z$", 'Interpreter', 'latex');
colororder("reef");
% ylim([-1 1]*maxError(3))
% subplot
subplot(3, 2, 2);
% plot(t, StateErrors(4, :));
scatter(t, StateErrors_CKF(4, :), 'x');
hold on
scatter(t, StateErrors_batch(4, :), '.');
Pplus2 = plot(t, 3*sig_4_CKF, '--');
Pminus2 = plot(t, -3*sig_4_CKF, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 3*sig_4_batch, '--');
Pminus2 = plot(t, -3*sig_4_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{x}$", 'Interpreter', 'latex');
colororder("reef");
legend("SRIF", "CKF", "SRIF 3$\sigma$", '', "CKF 3$\sigma$", '','Location', 'northeast', "Interpreter", "latex")
% ylim([-1 1]*maxError(4))
% subplot
subplot(3, 2, 4);
% plot(t, StateErrors(5, :));
scatter(t, StateErrors_CKF(5, :), 'x');
hold on
scatter(t, StateErrors_batch(5, :), '.');
Pplus1 = plot(t, 3*sig_5_CKF, '--');
Pminus1 = plot(t, -3*sig_5_CKF, '--');
Pminus1.SeriesIndex = Pplus1.SeriesIndex;
Pplus2 = plot(t, 3*sig_5_batch, '--');
Pminus2 = plot(t, -3*sig_5_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{y}$", 'Interpreter', 'latex');
colororder("reef");
% ylim([-1 1]*maxError(5))
% subplot
subplot(3, 2, 6);
% plot(t, StateErrors(6, :));
scatter(t, StateErrors_CKF(6, :), 'x');
hold on
scatter(t, StateErrors_batch(6, :), '.');
Pplus2 = plot(t, 3*sig_6_CKF, '--');
Pminus2 = plot(t, -3*sig_6_CKF, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
Pplus2 = plot(t, 3*sig_6_batch, '--');
Pminus2 = plot(t, -3*sig_6_batch, '--');
Pminus2.SeriesIndex = Pplus2.SeriesIndex;
xlabel('time [hours]');
ylabel("$\delta \dot{z}$", 'Interpreter', 'latex');
colororder("reef");
% ylim([-1 1]*maxError(6))
end

