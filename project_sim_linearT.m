clear all;
close all;
clc;
%% 
% Init
Vm = 10;
k1 = 5.5; k2 = 0.5;
Xm0 = -100; Ym0 = 20; Xt0 = 0; Yt0 = 0;
rho0 = 100;

psi_m0 = deg2rad(-45);
sigma0 = deg2rad(78.69);
theta0 = sigma0 + psi_m0;
psi_t = deg2rad(45);

% Vg = 1.5;
Vm_x0 = Vm * cos(psi_m0);
Vm_y0 = Vm * sin(psi_m0);
Wx = 0; Wy = 0; % wind velocity components

%%
tspan = [0 250];
initial_condition = [rho0, psi_m0, theta0, Xm0, Ym0,0,0];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, y] = ode45(@(t,y) tracksim(t, y, Vm,psi_t, Wx, Wy, k1, k2), tspan, initial_condition, options);

%%

% Calculate Vt for each time step after ODE solution
vt_sol = arrayfun(@(t) calculateVt(t), t);


rho = y(:,1);
psi_m = y(:,2);
theta = y(:,3);
Xm = y(:,4);
Ym = y(:,5);
Xt = y(:,6);
Yt = y(:,7);

%%
figure(1);
plot(t, vt_sol, 'B', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Vt (m/s)');
title('Vt over Time');
grid on;

figure(2);
plot(t, rho,'k','LineWidth',1.5);
title('rho vs t')
title('LOS distance vs time')
xlabel('t (s)');
ylabel('rho (m)');
grid on;

figure(3);
plot(Xt, Yt, 'r--');
hold on;
plot(Xm, Ym, 'k', 'LineWidth', 1.5);
title('Trajectory');
axis('equal');
xlabel(" X(m)");
ylabel("Y(m)");
grid on;

figure(4);
plot(t, rad2deg(psi_m),'k','LineWidth',1.5);
xlabel('t (s)');
ylabel('psi_m (degree)');
title('UAV heading angle vs time');
grid on;

% aN = (k1*theta)/(cosh(theta)-k2);
% figure(5);
% plot(t, aN,'k','LineWidth',1.5);
% title('Lateral Accelaration vs t')
% xlabel('t (s)');
% ylabel('aN');
% grid on;
%% Dynamics function for UAV
function dy = tracksim(t,y,Vm,psi_t,Wx, Wy,k1,k2)
    
    rho = y(1);
    psi_m = y(2);
    theta = y(3);
    Vt = calculateVt(t);


    Vg = sqrt(((Vm*cos(psi_m) + Wx)^2) + ((Vm * sin(psi_m) + Wy)^2));

    rho_dot = (Vt * cos(theta + psi_m - psi_t)) - Vg * cos(theta);


    aN = (k1*theta)/(cosh(theta)-k2);
    w = aN / Vg;
    psi_mdot = w;

    theta_dot = (((-Vt * sin(theta + psi_m - psi_t)) + (Vg * sin(theta)))/rho)-(psi_mdot);

    %state equation for trajectories
    Xm_dot = Vm * cos(psi_m) ;
    Ym_dot = Vm * sin(psi_m) ;

    Xt_dot = Vt * cos(psi_t);
    Yt_dot = Vt * sin(psi_t);
    
    dy = [rho_dot;psi_mdot;theta_dot;Xm_dot;Ym_dot;Xt_dot;Yt_dot];
end 
%%
function Vt = calculateVt(t)
    if t == 0
        Vt = 0;
    elseif t > 0 && t < 50
        Vt = 0.1 * t; % Increase linearly from 0 to 5
    elseif t >= 50 && t < 75
        Vt = 5; % Constant at 5
    elseif t >= 75 && t < 100
        % Linearly increase from 5 to 8
        slope = (8 - 5) / (100 - 75); % Calculate the slope
        Vt = 5 + slope * (t - 75); % Apply linear equation (y = mx + b)
    elseif t >= 100 && t < 125
        Vt = 8; % Constant at 8
    elseif t >= 125 && t < 200
        % Linearly decrease to 0 by time 200
        slope = -8 / (200 - 125); % Calculate the slope of the decrease
        Vt = 8 + slope * (t - 125); % Apply the linear equation
    elseif t >= 200
        Vt = 0; % Vt is 0 at time 200 and onwards
    else
        Vt = NaN; % Handles any unexpected values of t
    end
end

