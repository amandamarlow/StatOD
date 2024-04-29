function [] = plotDeltaX(t, deltaX, title)



% SC State
% figure
% subplot(3, 6, 1);
% sgtitle(title, 'Interpreter', 'latex');
% scatter(t, deltaX(1, :), '.')
% hold on
% xlabel('time [s]');
% ylabel("$\Delta x$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 6, 7);
% scatter(t, deltaX(2, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta y$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 6, 13);
% scatter(t, deltaX(3, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta z$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 6, 2);
% scatter(t, deltaX(4, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta \dot{x}$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 6, 8);
% scatter(t, deltaX(5, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta \dot{y}$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 6, 14);
% scatter(t, deltaX(6, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta \dot{z}$", 'Interpreter', 'latex');
% grid on

figure
% Constants
subplot(3, 2, 1);
scatter(t, deltaX(7, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta C_R$", 'Interpreter', 'latex');
grid on
%subplot
subplot(3, 2, 2);
scatter(t, deltaX(8, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta \mu E$", 'Interpreter', 'latex');
grid on
%subplot
subplot(3, 2, 3);
scatter(t, deltaX(9, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta \mu S$", 'Interpreter', 'latex');
grid on
% % subplot
% subplot(3, 2, 4);
% scatter(t, deltaX(10, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta P_{srAU}$", 'Interpreter', 'latex');
% grid on
% % subplot
% subplot(3, 2, 5);
% scatter(t, deltaX(11, :), '.');
% hold on
% xlabel('time [s]');
% ylabel("$\Delta A/m$", 'Interpreter', 'latex');
% grid on

% Station Locations
figure
subplot(3, 3, 1);
scatter(t, deltaX(end-8, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 4);
scatter(t, deltaX(end-7, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 7);
scatter(t, deltaX(end-6, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 2);
scatter(t, deltaX(end-5, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 5);
scatter(t, deltaX(end-4, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 8);
scatter(t, deltaX(end-3, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 3);
scatter(t, deltaX(end-2, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s3}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 6);
scatter(t, deltaX(end-1, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s3}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 3, 9);
scatter(t, deltaX(end, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s3}$", 'Interpreter', 'latex');
grid on

end

