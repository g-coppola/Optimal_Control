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

% Projection Matrices
B1 = [1 0 0 0; 0 1 0 0]; 
B2 = [0 0 1 0; 0 0 0 1];

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

color = [0, 0.447, 0.741]; % unico colore

set(gcf, 'Color', 'w');

titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

for i = 1:4
    subplot(4,1,i);
    plot(t, xs(:,i), 'Color', color, 'LineWidth', 2);
    
    title(titles{i}, 'Interpreter', 'latex', 'FontSize', 12)
    
    grid on;
    set(gca, ...
        'FontSize', 11, ...
        'FontName', 'Times New Roman', ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.2, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on', ...
        'LineWidth', 1);
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 13);
    end
end

sgtitle('\textbf{System Free Response}', 'Interpreter', 'latex', 'FontSize', 15);

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

colors = [0, 0.447, 0.741]; 
input_color = [0.1, 0.6, 0.2];

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

for i = 1:4
    nexttile((i-1)*2 + 1);
    
    plot(t, xs(i,:), 'Color', colors, 'LineWidth', 2);
    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
end

ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

plot(t, us, 'Color', input_color, 'LineWidth', 2);

title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

title(tlo, '\textbf{LQR (Infinite) - Zero Regulation}', ...
    'Interpreter', 'latex', 'FontSize', 15);


input("")
clear xs us t

% Predictive Control for Finite Horizon
[t, xs,us] = OC.OLQR(ft,x0,N);

% Plot
fig = figure('Units', 'inches'); 

colors = [0, 0.447, 0.741]; 
input_color = [0.1, 0.6, 0.2];

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

for i = 1:4
    nexttile((i-1)*2 + 1);
    
    plot(t, xs(i,:), 'Color', colors, 'LineWidth', 2);
    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
end

ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

plot(t, us, 'Color', input_color, 'LineWidth', 2);

title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

title(tlo, '\textbf{LQR (Model Predictive Finite) - Zero Regulation}', ...
    'Interpreter', 'latex', 'FontSize', 15);

input('Press to Polytopic Description (Unconstrained)')
close all
clear xs us t

%% POLYTOPIC DESCRIPTION (Unconstrained)
Nvals = 30;
xeq_vals = linspace(-0.4,0.4,Nvals);
param.tspan = [0 6];
N = round(param.tspan(2)/param.dt);

param.Q1 = diag([1 0.1 1 0.1]);
param.RR = 0.01*eye(param.m);
h = 1;

x_store = {};
t_store = {};
u_store = {};

[A_cell,B_cell,xeq_used] = build_pol(param,f,xeq_vals,x,u);
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

e = ell(QQ);

% Plot
fig = figure('Units', 'inches'); 

colors = [0, 0.447, 0.741]; 
input_color = [0.1, 0.6, 0.2];

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

for i = 1:4
    nexttile((i-1)*2 + 1);
    
    plot(t, xs(i,:), 'Color', colors, 'LineWidth', 2);
    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
end

ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

plot(t(1:end-1), us, 'Color', input_color, 'LineWidth', 2);

title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

title(tlo, '\textbf{Polytopic Description - Unconstrained}', ...
    'Interpreter', 'latex', 'FontSize', 15);

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
e.plotEllBoundary(B1,'k--');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(xs(1,:), xs(2,:), 'b', 'LineWidth', 1.8);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Unconstrained}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;

nexttile;
hold on;
e.plotEllBoundary(B2,'k--');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
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

x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;
input('Press to Polytopic Description (Constrained Input)')
close all
clear F t us xs Qe Q1 Q2

%% POLYTOPIC DESCRIPTION (Constrained Input)
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

ec = ell(QQ);

% Plot
fig = figure('Units', 'inches'); 

colors = [0, 0.447, 0.741]; 
input_color = [0.1, 0.6, 0.2];

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

for i = 1:4
    nexttile((i-1)*2 + 1);
    
    plot(t, xs(i,:), 'Color', colors, 'LineWidth', 2);
    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
end

ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

plot(t(1:end-1), us, 'Color', input_color, 'LineWidth', 2);
hold on
yline(u_max, '--k', 'LineWidth', 1.2);
yline(-u_max, '--k', 'LineWidth', 1.2);
legend({'$u(t)$ - Applied Torque [rad/s$^2$]','$u_{max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;
hold off;



title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-4.5 4.5])

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

title(tlo, '\textbf{Polytopic Description - Constrained Input}', ...
    'Interpreter', 'latex', 'FontSize', 15);

fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
ec.plotEllBoundary(B1,'--k');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(xs(1,:), xs(2,:), 'b', 'LineWidth', 1.8);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;

nexttile;
hold on;
ec.plotEllBoundary(B2,'--k');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
plot(param.x_eq(3), param.x_eq(4), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(xs(3,:), xs(4,:), 'b', 'LineWidth', 1.8);
xlabel('$x_3(t)$ - [rad]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_4(t)$ - [rad/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;

x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;

input('Comparison - Unconstrained Constrained')
close all
clear F t us xs Qe Q1 Q2

%% COMPARISON (Unconstrained vs Constrained)

fig = figure('Units', 'inches'); 

color_free   = [0, 0.447, 0.741];     % stati NON vincolati (blu)
color_constr = [0.850, 0.325, 0.098]; % stati vincolati (rosso)

input_free   = [0.1, 0.6, 0.2];       % input NON vincolato (verde)
input_constr = [0.494, 0.184, 0.556]; % input vincolato (viola)

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

% ===== STATI =====
for i = 1:4
    nexttile((i-1)*2 + 1);
    hold on;
    
    plot(t_store{1}, x_store{1}(i,:), ...
        '-', 'Color', color_free, 'LineWidth', 2);
    
    plot(t_store{2}, x_store{2}(i,:), ...
        '--', 'Color', color_constr, 'LineWidth', 2);
    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
    
    if i == 1
        legend({'Unconstrained', 'Constrained'}, ...
            'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best');
        legend boxoff;
    end
    
    hold off;
end

% ===== INPUT =====
ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

hold on;

% NON vincolato
plot(t_store{1}(1:end-1), u_store{1}, ...
    '-', 'Color', input_free, 'LineWidth', 2);

% vincolato
plot(t_store{2}(1:end-1), u_store{2}, ...
    '--', 'Color', input_constr, 'LineWidth', 2);

% limiti
yline(u_max, '--k', 'LineWidth', 1.2);
yline(-u_max, '--k', 'LineWidth', 1.2);

legend({'Unconstrained', 'Constrained', '$\pm u_{max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;

title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);

ylim([-4.5 4.5]);

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

hold off;

% ===== TITOLO =====
title(tlo, '\textbf{Polytopic Description - Constrained vs Unconstrained}', ...
    'Interpreter', 'latex', 'FontSize', 15);


fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');


nexttile;
hold on;
plot(x_store{1}(1,:), x_store{1}(2,:), 'b', 'LineWidth', 1.8);
plot(x_store{2}(1,:), x_store{2}(2,:), 'r', 'LineWidth', 1.8);
e.plotEllBoundary(B1,'--k')
ec.plotEllBoundary(B1,'--r');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
legend({'Unconstrained', 'Constrained'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;


nexttile;
hold on;
plot(x_store{1}(3,:), x_store{1}(4,:), 'b', 'LineWidth', 1.8);
plot(x_store{2}(3,:), x_store{2}(4,:), 'r', 'LineWidth', 1.8);
e.plotEllBoundary(B2,'--k')
ec.plotEllBoundary(B2,'--r');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('$x_3(t)$ - [rad]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_4(t)$ - [rad/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
legend({'Unconstrained', 'Constrained'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
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

ecw = ell(QQ);

% Plot
fig = figure('Units', 'inches'); 

color_free   = [0, 0.447, 0.741];     % stati NON vincolati (blu)
color_constr = [0.850, 0.325, 0.098]; % stati vincolati (rosso)

input_free   = [0.1, 0.6, 0.2];       % input NON vincolato (verde)
input_constr = [0.494, 0.184, 0.556]; % input vincolato (viola)

tlo = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

state_titles = {...
    '$x_1(t)$ - Ball Position [m]', ...
    '$x_2(t)$ - Ball Velocity [m/s]', ...
    '$x_3(t)$ - Beam Angle [rad]', ...
    '$x_4(t)$ - Beam Velocity [rad/s]'};

% ===== STATI =====
for i = 1:4
    nexttile((i-1)*2 + 1);
    hold on;
    
    plot(t, x_store{1}(i,:), ...
        '-', 'Color', color_free, 'LineWidth', 2);
    
    plot(t, xs(i,:), ...
        '--', 'Color', color_constr, 'LineWidth', 2);

    yline([-y_max(i), y_max(i)], '--k', 'LineWidth', 1.1, 'Alpha', 0.7);
    ylim([-y_max(i)-0.05 y_max(i)+0.05]); grid on;

    
    title(state_titles{i}, 'Interpreter', 'latex', 'FontSize', 11);
    
    grid on;
    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on', ...
        'GridAlpha', 0.15, ...
        'TickLabelInterpreter', 'latex', ...
        'Box', 'on');
    
    if i ~= 4
        set(gca, 'XTickLabel', []);
    else
        xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
    end
    
    if i == 1
        legend({'Unconstrained', 'Constrained','$\pm y_{max}$'}, ...
            'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best');
        legend boxoff;
    end
    
    hold off;
end

% ===== INPUT =====
ax_input = nexttile(2);
ax_input.Layout.TileSpan = [4 1];

hold on;

% NON vincolato
plot(t_store{1}(1:end-1), u_store{1}, ...
    '-', 'Color', input_free, 'LineWidth', 2);

% vincolato
plot(t(1:end-1), us, ...
    '--', 'Color', input_constr, 'LineWidth', 2);

% limiti
yline(u_max, '--k', 'LineWidth', 1.2);
yline(-u_max, '--k', 'LineWidth', 1.2);

legend({'Unconstrained', 'Constrained', '$\pm u_{max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;

title('$u(t)$ - Applied Torque', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u(t)$', 'Interpreter', 'latex', 'FontSize', 12);

ylim([-4.5 4.5]);

grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'GridAlpha', 0.15, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');

hold off;

% ===== TITOLO =====
title(tlo, '\textbf{Polytopic Description - Constrained vs Unconstrained}', ...
    'Interpreter', 'latex', 'FontSize', 15);



fig1 = figure('Units', 'inches');
tlo1 = tiledlayout(2, 1, 'TileSpacing', 'loose', 'Padding', 'compact');

nexttile;
hold on;
plot(x_store{1}(1,:), x_store{1}(2,:), 'b', 'LineWidth', 1.8);
plot(xs(1,:), xs(2,:), 'r', 'LineWidth', 1.8);
e.plotEllBoundary(B1,'--k')
ecw.plotEllBoundary(B1,'--r');
rectangle('Position',[-y_max(1), -y_max(2), 2*y_max(1), 2*y_max(2)], ...
          'EdgeColor', [0.85, 0.325, 0.098], 'LineStyle', '--', 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('$x_1(t)$ - [m]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_2(t)$ - [m/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_1, x_2)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
legend({'Unconstrained', 'Constrained'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;


nexttile;
hold on;
plot(x_store{1}(3,:), x_store{1}(4,:), 'b', 'LineWidth', 1.8);
plot(xs(3,:), xs(4,:), 'r', 'LineWidth', 1.8);
e.plotEllBoundary(B2,'--k')
ecw.plotEllBoundary(B2,'--r');
rectangle('Position',[-y_max(3), -y_max(4), 2*y_max(3), 2*y_max(4)], ...
          'EdgeColor', [0.85, 0.325, 0.098], 'LineStyle', '--', 'LineWidth', 1.2);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('$x_3(t)$ - [rad]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('$x_4(t)$ - [rad/s]', 'Interpreter', 'latex', 'FontSize', 13);
title('\boldmath{$(x_3, x_4)$} \textbf{Evolution - Constrained Input}', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on');
hold off;
legend({'Unconstrained', 'Constrained'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
legend boxoff;


