function [] = scatterDeltaXsc(t, deltaX, title)

figure

% SC State

subplot(3, 2, 1);
sgtitle(title, 'Interpreter', 'latex');
scatter(t, deltaX(1, :), '.')
hold on
xlabel('time [s]');
ylabel("$\Delta x$ [km]", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 2, 3);
scatter(t, deltaX(2, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta y$ [km]", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 2, 5);
scatter(t, deltaX(3, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta z$ [km]", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 2, 2);
scatter(t, deltaX(4, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{x}$ [km/s]", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 2, 4);
scatter(t, deltaX(5, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{y}$ [km/s]", 'Interpreter', 'latex');
grid on
% subplot
subplot(3, 2, 6);
scatter(t, deltaX(6, :), '.');
hold on
xlabel('time [s]');
ylabel("$\Delta \dot{z}$ [km/s]", 'Interpreter', 'latex');
grid on

end

