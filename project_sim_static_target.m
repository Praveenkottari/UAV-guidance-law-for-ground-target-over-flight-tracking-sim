clear all;
close all;
clc;
%% 
%Init
Vm = 10;
Vt =0;

k1 = 5.5; k2 = 0.5;

Xm0 = -100; Ym0 = 20; Xt0 = 0; Yt0 = 0;
rho0 = 101.98;

psi_m0 = deg2rad(-56.31);
psi_t = deg2rad(30);
sigma0 = deg2rad((-101.31));
theta0 = deg2rad((45));

Wx = 0; Wy = 0; % wind velocity components

%%
tspan = [0 200];
initial_condition = [rho0,psi_m0,theta0,Xm0,Ym0];

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t,y] = ode45(@(t,y) tracksim(y,Vm,Vt,psi_t,Wx, Wy,k1,k2), tspan, initial_condition, options);

%%
rho = y(:,1);
psi_m = y(:,2);
theta = y(:,3);
Xm = y(:,4);
Ym = y(:,5);
%% Plots
figure(1);
plot(t, rho,'k','LineWidth',1.5);
xlabel('t (s)');
ylabel('rho (m)');
title('LOS distance vs time');
grid on;

%trajectory
Xt = Xt0 + zeros(length(t),1);
Yt = Yt0 + zeros(length(t),1);
figure(2);
plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
hold on;
plot(Xm,Ym,'k','LineWidth',1.5);
title('Trajectory');
axis('equal');
xlabel(" X(m)")
ylabel("Y(m)")
%axis([-100 80 -80 80])
grid on;

figure(3);
plot(t, rad2deg(psi_m),'k','LineWidth',1.5);
xlabel('t (s)');
ylabel('psi_m (degree)');
title('UAV heading angle vs time');
grid on;

% aN = (k1*theta)/(cosh(theta)-k2);
% figure(4);
% plot(t, aN,'k','LineWidth',1.5);
% title('Lateral Accelaration vs t')
% xlabel('t (s)');
% ylabel('aN');
% grid on;
%%
function dy = tracksim(y,Vm,Vt,psi_t,Wx, Wy,k1,k2)
    
    rho = y(1);
    psi_m = y(2);
    theta = y(3);

    Vg = sqrt((((Vm*cos(psi_m)) + Wx)^2) + (((Vm * sin(psi_m)) + Wy)^2));

    rho_dot = (Vt * cos(theta + psi_m - psi_t)) - Vg * cos(theta);
   
    aN = (k1*theta)/(cosh(theta)-k2);
    w = aN / Vg;
    psi_mdot = w;

    theta_dot = (((-Vt * sin(theta + psi_m - psi_t)) + (Vg * sin(theta)))/rho)-(psi_mdot);

    %state equation for trajectories
    Xm_dot = Vm * cos(psi_m) + Wx;
    Ym_dot = Vm * sin(psi_m) + Wy ;

    dy = [rho_dot;psi_mdot;theta_dot;Xm_dot;Ym_dot];
end 