function [] = scatterDeltaX(t, deltaX, title)

figure

% SC State

subplot(3, 6, 1);
sgtitle(title, 'Interpreter', 'latex');
scatter(t, deltaX(1, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 7);
scatter(t, deltaX(2, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 13);
scatter(t, deltaX(3, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 2);
scatter(t, deltaX(4, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{x}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 8);
scatter(t, deltaX(5, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{y}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 14);
scatter(t, deltaX(6, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{z}$", 'Interpreter', 'latex');
grid on

% Constants

subplot(3, 6, 3);
scatter(t, deltaX(7, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta \mu$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 9);
scatter(t, deltaX(8, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta J2$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 15);
scatter(t, deltaX(9, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta C_D$", 'Interpreter', 'latex');
grid on

% Station Locations

subplot(3, 6, 4);
scatter(t, deltaX(10, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 10);
scatter(t, deltaX(11, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 16);
scatter(t, deltaX(12, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s1}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 5);
scatter(t, deltaX(13, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 11);
scatter(t, deltaX(14, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 17);
scatter(t, deltaX(15, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s2}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 6);
scatter(t, deltaX(16, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x_{s3}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 12);
scatter(t, deltaX(17, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y_{s3}$", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 6, 18);
scatter(t, deltaX(18, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z_{s3}$", 'Interpreter', 'latex');
grid on

end

