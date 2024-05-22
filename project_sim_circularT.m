clear all;
close all;
clc;
%% 
%Init
Vm = 10;
Vt =5;
aT = 0.05;

k1 = 5.5; k2 = 0.5;

Xm0 = -100; Ym0 = 20; Xt0 = 0; Yt0 = 0;

rho0 = 101.98;

psi_m0 = deg2rad(-45);

theta0 = deg2rad(33.69);

Wx = 0; Wy = 3; % wind velocity components

%%
tspan = [0 550];
initial_condition = [rho0,psi_m0,theta0,Xm0,Ym0,0.01,0,0];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t,y] = ode45(@(t,y) tracksim(y,Vm,Vt,aT,Wx, Wy,k1,k2), tspan, initial_condition, options);

%%
rho = y(:,1);
psi_m = y(:,2);
theta = y(:,3);
Xm = y(:,4);
Ym = y(:,5);
Xt = y(:,7);
Yt = y(:,8);

%% Plots

figure(1);
plot(t, rho,'k','LineWidth',1.5);
title('rho vs t')
xlabel('t (s)');
ylabel('rho (m)');
grid on;

figure(2);
plot(Xt, Yt, 'r--','LineWidth',1.5);
hold on;
plot(Xm,Ym,'k','LineWidth',1.5);
title('Trajectory');
axis('equal');
xlabel(" X(m)")
ylabel("Y(m)")
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
function dy = tracksim(y,Vm,Vt,aT,Wx, Wy,k1,k2)
    
    rho = y(1);
    psi_m = y(2);
    theta = y(3);
    psi_t = y(6);
    

    Vg = sqrt(((Vm*cos(psi_m) + Wx)^2) + ((Vm * sin(psi_m) + Wy)^2));

    rho_dot = (Vt * cos(theta + psi_m - psi_t)) - Vg * cos(theta);
   
    aN = (k1*theta)/(cosh(theta)-k2);
    w = aN / Vg;
    psi_mdot = w;

    theta_dot = (((-Vt * sin(theta + psi_m - psi_t)) + (Vg * sin(theta)))/rho)-(psi_mdot);

    %state equation for trajectories
    Xm_dot = Vg * cos(psi_m) + Wx;
    Ym_dot = Vg * sin(psi_m) + Wy ;
    Xt_dot = Vt * cos(psi_t) ;
    Yt_dot = Vt * sin(psi_t)  ;


    psi_tdot = aT/Vt;
    
    dy = [rho_dot;psi_mdot;theta_dot;Xm_dot;Ym_dot;psi_tdot;Xt_dot;Yt_dot];
end 