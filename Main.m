%% BALL AND BEAM    
clear
close all
clc

%% DYNAMICAL SYSTEM

% Parameters
param.I = 0.5;
param.g = 9.81;
param.M = 2;
param.J = 25e-3;
param.r = 0.05;

% State/Input Dimension
param.m = 1;
param.n = 4;

% Sampling Time
param.dt = 0.1;

param.tspan = [0 8];

% System
syms u;                
syms x [4 1]; 

I = param.I;
g = param.g;
m = param.M;
J = param.J;
r = param.r;
param.MM = m+(J/r^2);
M = param.MM;

f = [x(2);
    (1/M)*(m*x(1)*x(4)^2-m*g*sin(x(3)));
    x(4);
    (1/(I+(m*x(1)^2)))*(u-2*m*x(1)*x(2)*x(4)-m*g*x(1)*cos(x(3)))
    ];

% Equilibrium
param.x_eq = [0 0 0 0]';
param.u_eq = 0;

% Initial Conditions
x0 = [0.15 0 0 0]';

%% FREE RESPONSE
f_fun = matlabFunction(f,'Vars',{x,u});

[t,xs] = ode45(@(t,x)f_fun(x,0),param.tspan,x0);

fig = figure('Units', 'inches');
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.466, 0.674, 0.188; 0.494, 0.184, 0.556];

hold on;
set(gca, 'ColorOrder', colors, 'FontSize', 11, 'FontName', 'Times New Roman');
plot(t, xs, 'LineWidth', 2);

xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('State variables $x(t)$', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{System Free Response}', 'Interpreter', 'latex', 'FontSize', 15);

lgd = legend({'$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', ...
              '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]'}, ...
             'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;

grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.2, ...
    'TickLabelInterpreter', 'latex', 'Box', 'on', 'LineWidth', 1);
hold off;

input('Press to LQR (Zero Regulation)')
close all   
clear xs t

%% LQR (Zero Regulation)
N = round(param.tspan(2)/param.dt);

param.Q = diag([10 1 10 1]);
param.R = 0.01*eye(param.m);
param.S = eye(param.n);

ft = @(t, x_val, u_val) f_fun(x_val, u_val);

OC = OptControl(param,f,x,u);

% Control for infinite Horizon
[t, xs, us, e] = OC.lqr(ft,x0,N);

% Plot
fig = figure('Units', 'inches'); 
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.466, 0.674, 0.188; 0.494, 0.184, 0.556];

tlo = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
set(gca, 'ColorOrder', colors, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, xs, 'LineWidth', 2);
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');


lgd1 = legend({'$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', ...
               '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]'}, ...
              'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;
hold off;

nexttile;
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, us, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');

lgd2 = legend({'$u(t)$ - Applied Torque [rad/s$^2$]'}, ...
              'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;
hold off;

title(tlo, '\textbf{LQR (Infinite) - Zero Regulation}', 'Interpreter', 'latex', 'FontSize', 15);

input("")
clear xs us t

% Predictive Control for Finite Horizon
[t, xs,us] = OC.OLQR(ft,x0,N);

% Plot
fig = figure('Units', 'inches');
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.466, 0.674, 0.188; 0.494, 0.184, 0.556];

tlo = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile;
hold on;
set(gca, 'ColorOrder', colors, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, xs, 'LineWidth', 2);
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', ...
        '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;
hold off;

nexttile;
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, us, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$u(t)$ - Applied Torque [rad/s$^2$]'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;
hold off;

title(tlo, '\textbf{LQR (Model Predictive Finite) - Zero Regulation}', 'Interpreter', 'latex', 'FontSize', 15);

input('Press to Polytopic Description (Unconstrained)')
close all
clear xs us t

%% POLYTOPIC DESCRIPTION (Unconstrained)
xeq_vals = {0 0.2 0.4};
param.tspan = [0 8];
N = round(param.tspan(2)/param.dt);

param.Q1 = diag([1 0.1 1 0.1]);
param.RR = 0.01*eye(param.m);
h = 1;

x_store = {};
t_store = {};
u_store = {};

[A_cell,B_cell] = build_pol(param,f,xeq_vals,x,u);
C = eye(param.n);

mpc = MPC(f_fun,ft,param,A_cell,B_cell,C);

xs = zeros(param.n,N);
us = zeros(param.m,N);
xs(:,1) = x0;


xk = x0;
for k = 1:N/h
    [F, Qe] = mpc.unconstrainedRHC(xk);
    if k == 1
        QQ = Qe;
    end
    for i = 1:h
        uk = F*xk;
        dx = ft(0,xk,uk);
        xk1 = xk+param.dt*dx;
        xs(:,k+i) = xk1;
        us(:,k+i-1) = uk;
        xk = xk1;
    end
end

t = 0:param.dt:N*param.dt;

Q1 = QQ(1:2,1:2);
Q2 = QQ(3:4,3:4);

e1 = ell(Q1);
e2 = ell(Q2);

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
e1.bestPlot([0.65, 0.65, 0.65]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(param.x_eq(1), param.x_eq(2), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(xs(1,:), xs(2,:), 'b', 'LineWidth', 1.8);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Unconstrained}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;

nexttile;
hold on;
e2.bestPlot([0.65, 0.65, 0.65]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(param.x_eq(3), param.x_eq(4), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(xs(3,:), xs(4,:), 'b', 'LineWidth', 1.8);
xlabel('$x_3(t)$ - [rad]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_4(t)$ - [rad/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Unconstrained}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
 
fig2 = figure('Units', 'inches');
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.466, 0.674, 0.188; 0.494, 0.184, 0.556];
tlo2 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile;
hold on;
set(gca, 'ColorOrder', colors, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, xs, 'LineWidth', 2);
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', ...
        '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;
hold off;

nexttile;
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t(1:end-1), us, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$u(t)$ - Applied Torque [rad/s$^2$]'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;
hold off;

title(tlo2, '\textbf{Polytopic Description - Unconstrained}', 'Interpreter', 'latex', 'FontSize', 15);

x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;
input('Press to Polytopic Description (Constrained Input)')
close all
clear F t us xs Qe Q1 Q2

% POLYTOPIC DESCRIPTION (Constrained Input)
u_max = 4 + param.u_eq;

xs = zeros(param.n,N);
us = zeros(param.m,N);
xs(:,1) = x0;


xk = x0;
for k = 1:N/h
    [F, Qe] = mpc.constrainedInputRHC(xk,u_max);
    if k == 1
        QQ = Qe;
    end
    for i = 1:h
        uk = F*xk;
        dx = ft(0,xk,uk);
        xk1 = xk+param.dt*dx;
        xs(:,k+i) = xk1;
        us(:,k+i-1) = uk;
        xk = xk1;
    end
end

t = 0:param.dt:N*param.dt;
Q1 = QQ(1:2,1:2);
Q2 = QQ(3:4,3:4);

e1c = ell(Q1);
e2c = ell(Q2);

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
e1.bestPlot([0.65, 0.65, 0.65]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(param.x_eq(1), param.x_eq(2), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(xs(1,:), xs(2,:), 'b', 'LineWidth', 1.8);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;

nexttile;
hold on;
e2.bestPlot([0.65, 0.65, 0.65]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(param.x_eq(3), param.x_eq(4), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(xs(3,:), xs(4,:), 'b', 'LineWidth', 1.8);
xlabel('$x_3(t)$ - [rad]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_4(t)$ - [rad/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
 
fig2 = figure('Units', 'inches');
colors = [0, 0.447, 0.741; 0.85, 0.325, 0.098; 0.466, 0.674, 0.188; 0.494, 0.184, 0.556];
tlo2 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile;
hold on;
set(gca, 'ColorOrder', colors, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t, xs, 'LineWidth', 2);
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', ...
        '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;
hold off;

nexttile;
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t(1:end-1), us, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);
yline(u_max, '--k', 'LineWidth', 1.2);
yline(-u_max, '--k', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridAlpha', 0.15, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'$u(t)$ - Applied Torque [rad/s$^2$]','$u_{max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;
hold off;

title(tlo2, '\textbf{Polytopic Description - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 15);

x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;

input('Comparison - Unconstrained Constrained')
close all
clear F t us xs Qe Q1 Q2

%% COMPARISON (Unconstrained vs Constrained)

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(1,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t_store{2}, x_store{2}(1,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
ylabel('Position $x_1$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(2,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t_store{2}, x_store{2}(2,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
ylabel('Velocity $x_2$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');

nexttile([1 2]); hold on;
e1.bestPlot([0.65 0.65 0.65]); 
e1c.bestPlot([0.9 0.4 0.4]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
p1 = plot(x_store{1}(1,:), x_store{1}(2,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 2);
p2 = plot(x_store{2}(1,:), x_store{2}(2,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 2);
plot(param.x_eq(1), param.x_eq(2), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
xlabel('$x_1$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend([p1 p2], {'Unconstrained', 'Constrained'}, 'Interpreter', 'latex', 'Location', 'northeast');
legend boxoff;

title(tlo1, '\textbf{Ball States Evolution: Unconstrained vs Constrained}', 'Interpreter', 'latex', 'FontSize', 14);


fig2 = figure('Units', 'inches');
tlo2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(3,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t_store{2}, x_store{2}(3,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
ylabel('Angle $x_3$ [rad]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(4,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t_store{2}, x_store{2}(4,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
ylabel('Ang. Vel. $x_4$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');

nexttile([1 2]); hold on;
e2.bestPlot([0.65 0.65 0.65]);
e2c.bestPlot([0.9 0.4 0.4]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
p3 = plot(x_store{1}(3,:), x_store{1}(4,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 2);
p4 = plot(x_store{2}(3,:), x_store{2}(4,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 2);
plot(param.x_eq(3), param.x_eq(4), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
xlabel('$x_3$ [rad]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_4$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
grid on; set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend([p3 p4], {'Unconstrained', 'Constrained'}, 'Interpreter', 'latex', 'Location', 'northeast');
legend boxoff;

title(tlo2, '\textbf{Beam States Evolution: Unconstrained vs Constrained}', 'Interpreter', 'latex', 'FontSize', 14);

fig3 = figure('Units', 'inches');
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(t_store{1}(1:end-1), u_store{1}, 'Color', [0, 0.447, 0.741], 'LineWidth', 2);
plot(t_store{2}(1:end-1), u_store{2}, 'Color', [0.85, 0.325, 0.098], 'LineWidth', 2);
yline(u_max, '--k', 'LineWidth', 1.2, 'Alpha', 0.6);
yline(-u_max, '--k', 'LineWidth', 1.2, 'Alpha', 0.6);

xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Control Input Comparison}', 'Interpreter', 'latex', 'FontSize', 15);
ylim([-5 5]);
grid on; set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'TickLabelInterpreter', 'latex', 'Box', 'on');

legend({'Unconstrained', 'Constrained', 'Saturation'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
legend boxoff;

input('Component Wise Constraints')
close all

%% POLYTOPIC DESCRIPTION (Component Wise Constraints)
y_max = [0.4, 0.25, 0.6, 1]';

xs = zeros(param.n,N);
us = zeros(param.m,N);
xs(:,1) = x0;


xk = x0;
for k = 1:N/h
    [F, Qe] = mpc.compWiseOutputRHC(xk,u_max,y_max);
    if k == 1
        QQ = Qe;
    end
    for i = 1:h
        uk = F*xk;
        dx = ft(0,xk,uk);
        xk1 = xk+param.dt*dx;
        xs(:,k+i) = xk1;
        us(:,k+i-1) = uk;
        xk = xk1;
    end
end

t = 0:param.dt:N*param.dt;

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(1,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t, xs(1,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
yline([-y_max(1), y_max(1)], '--k', 'LineWidth', 1.1, 'Alpha', 0.7);
ylim([-0.5 0.5]); grid on;
ylabel('Position $x_1$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'Unconstrained', 'Constrained', 'Limits'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
legend boxoff;

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(2,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t, xs(2,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
yline([-y_max(2), y_max(2)], '--k', 'LineWidth', 1.1, 'Alpha', 0.7);
ylim([-0.5 0.5]); grid on;
ylabel('Velocity $x_2$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'Unconstrained', 'Constrained', 'Limits'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
legend boxoff;

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(3,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t, xs(3,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
yline([-y_max(3), y_max(3)], '--k', 'LineWidth', 1.1, 'Alpha', 0.7);
ylim([-0.7 0.7]); grid on;
ylabel('Angle $x_3$ [rad]', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'Unconstrained', 'Constrained', 'Limits'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
legend boxoff;

nexttile; hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
plot(t_store{1}, x_store{1}(4,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(t, xs(4,:), 'Color', [0.85, 0.325, 0.098], 'LineWidth', 1.8);
yline([-y_max(4), y_max(4)], '--k', 'LineWidth', 1.1, 'Alpha', 0.7);
ylim([-1.1 1.1]); grid on;
ylabel('Ang. Vel. $x_4$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
legend({'Unconstrained', 'Constrained', 'Limits'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
legend boxoff;

title(tlo1, '\textbf{State Evolution under State and Input Constraints}', 'Interpreter', 'latex', 'FontSize', 15);

fig2 = figure('Units', 'inches');
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);

plot(t(1:end-1), us, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);

yline([u_max, -u_max], '--k', 'LineWidth', 1.3, 'Alpha', 0.8);

xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('Input $u(t)$ [rad/s$^2$]', 'Interpreter', 'latex', 'FontSize', 13);
title('\textbf{Constrained Input Signal}', 'Interpreter', 'latex', 'FontSize', 15);
ylim([-4.5 4.5]);
grid on;
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'TickLabelInterpreter', 'latex', 'Box', 'on');

legend({'$u(t)$ - Applied Torque [rad/s$^2$]','$u_{max}$'}, 'Interpreter', 'latex', 'Location', 'best');
legend boxoff;


% ELLIPSOID (Component Wise Constraints)

Q1 = QQ(1:2,1:2);
Q2 = QQ(3:4,3:4);

e1cc = ell(Q1);
e2cc = ell(Q2);

fig = figure('Units', 'inches');
tlo = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile;
hold on;

e1.bestPlot([0.65, 0.65, 0.65], 1); 
e1cc.bestPlot([0.9, 0.4, 0.4]);     

rectangle('Position',[-y_max(1), -y_max(2), 2*y_max(1), 2*y_max(2)], ...
          'EdgeColor', [0.85, 0.325, 0.098], 'LineStyle', '--', 'LineWidth', 1.2);

p_un = plot(x_store{1}(1,:), x_store{1}(2,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
p_co = plot(xs(1,:), xs(2,:), 'Color', [0.850, 0.325, 0.098], 'LineWidth', 1.8);

plot(param.x_eq(1), param.x_eq(2), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'TickLabelInterpreter', 'latex', 'Box', 'on');
xlabel('$x_1(t)$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2(t)$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Constrained Input and State}', 'Interpreter', 'latex', 'FontSize', 13);
grid on;

legend([p_un, p_co], {'Unconstrained Trajectory', 'Constrained Trajectory'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9);
legend boxoff;
hold off;

nexttile;
hold on;
e2.bestPlot([0.65, 0.65, 0.65], 1);
e2cc.bestPlot([0.9, 0.4, 0.4]);

rectangle('Position',[-y_max(3), -y_max(4), 2*y_max(3), 2*y_max(4)], ...
          'EdgeColor', [0.85, 0.325, 0.098], 'LineStyle', '--', 'LineWidth', 1.2);

plot(x_store{1}(3,:), x_store{1}(4,:), 'Color', [0, 0.447, 0.741], 'LineWidth', 1.8);
plot(xs(3,:), xs(4,:), 'Color', [0.850, 0.325, 0.098], 'LineWidth', 1.8);
plot(param.x_eq(3), param.x_eq(4), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'TickLabelInterpreter', 'latex', 'Box', 'on');
xlabel('$x_3(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_4(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Constrained Input and State}', 'Interpreter', 'latex', 'FontSize', 13);
grid on;

legend({'Invariant Set (Unconst.)', 'Invariant Set (Const.)'}, ...
       'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9);
legend boxoff;
hold off;