clear
clc
close all

X0 = [1; 0]; % [m,m/s]
dt = 0.1;
tspan = 0:dt:6.2;
tm1 = 1.5;
targ1 = [-2, 2]; % [r_targ, t_targ]
tm2 = 2;
targ2 = [-2, 5]; % [r_targ, t_targ]
tm3 = 5;
targ3 = [0, 6.2]; % [r_targ, t_targ]
targs = [tm1, targ1; tm2, targ2; tm3, targ3];

sigDV = 0.1; % m/s

Dtrue0 = [0.1^2, 0.0099; 0.0099, 0.1^2]; % true dispersion
Pnav0 = [0.02^2, 0; 0, 0.01^2]; % nav

M = eye(2);
C0 = [Dtrue0, Dtrue0*M'; M*Dtrue0, M*Dtrue0*M'+Pnav0];

% % options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% S0 = [X0, reshape(eye(2), [], 1)];
% [~,S] = ode45(@(t,S) hw2ODE(t,S,n,constants), tspan, S0);% , options);

B = [zeros(1); eye(1)];
H = [1, 0];
% G will be 1by2
S = diag([0, 0.1^2]); % equivalent to mapping with B and M
% R = 0.01^2;
R = .1^2;
%% Reference Trajectory

Xref = [X0, zeros(2,length(tspan)-1)];
dVhist = zeros(1,3);
for k = 2:length(tspan)
    t = tspan(k);
    Phi = STM(tspan(k)-tspan(k-1));
    Xmin = Phi*Xref(:,k-1);
    if t == tm1
        targVec = targ1;
        dV = waypoint(Xmin, t, targVec);
        dVhist(1) = dV;
    elseif t == tm2
        targVec = targ2;
        dV = waypoint(Xmin, t, targVec);
        dVhist(2) = dV;
    elseif t == tm3
        targVec = targ3;
        dV = waypoint(Xmin, t, targVec);   
        dVhist(3) = dV;
    else
        dV = 0;
    end
    Xref(:,k) = Xmin + B*dV;

end

figure
subplot(2,1,1)
scatter(tspan, Xref(1,:))
ylabel("position")
subplot(2,1,2)
scatter(tspan, Xref(2,:))
xlabel("time")
ylabel("velocity")
%% Monte Carlo

N = 1000;
Xtrue0 = mvnrnd(X0, Dtrue0, N)';
Xnav0 = Xtrue0 + mvnrnd(zeros(size(X0)), Pnav0, N)';
[Xtrue, Xnav] = MonteCarlo(tspan, Xtrue0, Xnav0, Pnav0, targs, S, R, N);
%% LINCOV

[C, Pnav] = LINCOV(tspan, X0, Dtrue0, Pnav0, H, S, R, targs);
Dtrue = zeros(2,2,length(tspan));
Dnav = zeros(2,2,length(tspan));
PtrueNav = zeros(2,2,length(tspan));
for i = 1:length(tspan)
    Dtrue(:,:,i) = [eye(2), zeros(2)]*C(:,:,i)*[eye(2), zeros(2)]';
    Dnav(:,:,i) = [zeros(2), eye(2)]*C(:,:,i)*[zeros(2), eye(2)]';
    PtrueNav(:,:,i) = [eye(2), -eye(2)]*C(:,:,i)*[eye(2), -eye(2)]';
end

%% Plotting
grayColor = [.7 .7 .7];

figure
scatter(Xtrue0(1,:),Xtrue0(2,:), '.')
hold on
scatter(Xnav0(1,:),Xnav0(2,:), '.')
legend("true", "nav")
xlabel("Initial Position")
ylabel("Initial Velocity")
axis equal

sig_rNav = squeeze(Pnav(1,1,:).^(1/2))';
sig_vNav = squeeze(Pnav(2,2,:).^(1/2))';
sig_rNavTrue = squeeze(PtrueNav(1,1,:).^(1/2))';
sig_vNavTrue = squeeze(PtrueNav(2,2,:).^(1/2))';
sigD_rNav = squeeze(Dnav(1,1,:).^(1/2))';
sigD_vNav = squeeze(Dnav(2,2,:).^(1/2))';
sigD_rTrue = squeeze(Dtrue(1,1,:).^(1/2))';
sigD_vTrue = squeeze(Dtrue(2,2,:).^(1/2))';

% figure
% subplot(2,1,1)
% hold on
% for i = 1:N
%     plot(tspan, squeeze(Xnav(1,:,i)), 'Color', grayColor)
% end
% plot(tspan,Xref(1,:), '--')
% plot(tspan, Xref(1,:)+3*sigD_rTrue)
% plot(tspan, Xref(1,:)-3*sigD_rTrue)
% ylabel("position")
% subplot(2,1,2)
% hold on
% for i = 1:N
%     plot(tspan, squeeze(Xnav(2,:,i)), 'Color', grayColor)
% end
% plot(tspan,Xref(2,:), '--')
% plot(tspan, Xref(2,:)+3*sigD_vTrue)
% plot(tspan, Xref(2,:)-3*sigD_vTrue)
% xlabel("time")
% ylabel("velocity")
% sgtitle("Trajectory")

figure
subplot(2,1,1)
hold on
for i = 1:N
    plot(tspan, squeeze(Xtrue(1,:,i)), 'Color', grayColor)
end
plot(tspan,Xref(1,:), '--')
plot(tspan, Xref(1,:)+3*sigD_rTrue)
plot(tspan, Xref(1,:)-3*sigD_rTrue)
ylabel("position")
subplot(2,1,2)
hold on
for i = 1:N
    plot(tspan, squeeze(Xtrue(2,:,i)), 'Color', grayColor)
end
plot(tspan,Xref(2,:), '--')
plot(tspan, Xref(2,:)+3*sigD_vTrue)
plot(tspan, Xref(2,:)-3*sigD_vTrue)
xlabel("time")
ylabel("velocity")
sgtitle("True Dispersion")


figure
subplot(2,1,1)
hold on
for i = 1:N
    plot(tspan, squeeze(Xnav(1,:,i))-Xref(1,:), 'Color', grayColor)
end
plot(tspan, 3*sigD_rNav)
plot(tspan, -3*sigD_rNav)
ylabel("position")
subplot(2,1,2)
hold on
for i = 1:N
    plot(tspan, squeeze(Xnav(2,:,i))-Xref(2,:), 'Color', grayColor)
end
plot(tspan, 3*sigD_vNav)
plot(tspan, -3*sigD_vNav)
xlabel("time")
ylabel("velocity")
sgtitle("Nav Dispersion")

% figure
% subplot(2,1,1)
% hold on
% for i = 1:N
%     plot(tspan, squeeze(Xnav(1,:,i))-squeeze(Xtrue(1,:,i)), 'Color', grayColor)
% end
% plot(tspan, 3*sig_rNav)
% plot(tspan, -3*sig_rNav)
% ylabel("position")
% subplot(2,1,2)
% hold on
% for i = 1:N
%     plot(tspan, squeeze(Xnav(2,:,i))-squeeze(Xtrue(2,:,i)), 'Color', grayColor)
% end
% plot(tspan, 3*sig_vNav)
% plot(tspan, -3*sig_vNav)
% xlabel("time")
% ylabel("velocity")
% sgtitle("True Nav Error")
figure
subplot(2,1,1)
hold on
for i = 1:N
    plot(tspan, squeeze(Xnav(1,:,i))-squeeze(Xtrue(1,:,i)), 'Color', grayColor)
end
plot(tspan, 3*sig_rNavTrue)
plot(tspan, -3*sig_rNavTrue)
ylabel("position")
subplot(2,1,2)
hold on
for i = 1:N
    plot(tspan, squeeze(Xnav(2,:,i))-squeeze(Xtrue(2,:,i)), 'Color', grayColor)
end
plot(tspan, 3*sig_vNavTrue)
plot(tspan, -3*sig_vNavTrue)
xlabel("time")
ylabel("velocity")
sgtitle("True Nav Error")