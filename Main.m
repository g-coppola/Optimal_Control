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
param.MM = m+J/r^2;
M = param.MM;

f = [x(2);
    -(m*g/M)*sin(x(3))+(m/M)*x(1)*x(4)^2;
    x(4);
    (1/(m+x(1)^2))*(u-2*m*x(1)*x(2)*x(4)-m*g*x(1)*cos(x(3)))
    ];

% Equilibrium
param.x_eq = [0 0 0 0]';
param.u_eq = 0;

% Initial Conditions
x0 = param.x_eq + [0.15 0 0 0]';

%% FREE RESPONSE
f_fun = matlabFunction(f,'Vars',{x,u});

[t,xs] = ode45(@(t,x)f_fun(x,0),param.tspan,x0);

figure;
plot(t,xs,LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
grid
sgtitle('Free Response','Interpreter','latex')

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

figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t,us,'-g',LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]', 'Interpreter', 'latex', 'Location', 'best')
grid

sgtitle('LQR (Infinite) - Zero Regulation', 'Interpreter', 'latex')

input("")
clear xs us t

% Predictive Control for Finite Horizon
[t, xs,us] = OC.OLQR(ft,x0,N);

figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t,us,'-g',LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]', 'Interpreter', 'latex', 'Location', 'best')
grid

sgtitle('LQR (Model Predictive Finite) - Zero Regulation', 'Interpreter','latex')

input('Press to Polytopic Description (Unconstrained)')
close all
clear xs us t

%% POLYTOPIC DESCRIPTION (Unconstrained)
xeq_vals = {-0.4 0 0.4};
param.tspan = [0 7];
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

figure(1)
subplot(2,1,1)
e1.bestPlot([1 0 0])
hold on
plot(param.x_eq(1),param.x_eq(2),'ok','LineWidth',1)
plot(xs(1,:),xs(2,:),'b','LineWidth',1.5);
hold off
title('\big($x_1(t)$,$x_2(t)$\big) Evolution - Unconstrained','Interpreter','latex')
xlabel('$x_1(t)$ - [m]','Interpreter','latex')
ylabel('$x_2(t)$ - [m/s]','Interpreter','latex')
grid


subplot(2,1,2)
e2.bestPlot([1 0 0])
hold on
plot(param.x_eq(3),param.x_eq(4),'ok','LineWidth',1)
plot(xs(3,:),xs(4,:),'b','LineWidth',1.5);
hold off
title('\big($x_3(t)$,$x_4(t)$\big) Evolution - Unconstrained','Interpreter','latex')
xlabel('$x_3(t)$ - [rad]','Interpreter','latex')
ylabel('$x_3(t)$ - [rad/s]','Interpreter','latex')
grid
 
figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t(1:end-1),us,'-g',LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]', 'Interpreter', 'latex', 'Location', 'best')
grid

sgtitle('Polytopic Descrption - Unconstrained','Interpreter','latex')


x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;
input('Press to Polytopic Description (Constrained Input)')
close all
clear F t us xs Qe Q1 Q2

% POLYTOPIC DESCRIPTION (Constrained Input)
u_max = 3.5 + param.u_eq;

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

figure(1)
subplot(2,1,1)
e1c.bestPlot([1 0 0])
hold on
plot(param.x_eq(1),param.x_eq(2),'ok','LineWidth',1)
plot(xs(1,:),xs(2,:),'b','LineWidth',1.5);
hold off
title('\big($x_1(t)$,$x_2(t)$\big) Evolution - Constrained','Interpreter','latex')
xlabel('$x_1(t)$ - [m]','Interpreter','latex')
ylabel('$x_2(t)$ - [m/s]','Interpreter','latex')
grid

subplot(2,1,2)
e2c.bestPlot([1 0 0])
hold on
plot(param.x_eq(3),param.x_eq(4),'ok','LineWidth',1)
plot(xs(3,:),xs(4,:),'b','LineWidth',1.5);
hold off
title('\big($x_3(t)$,$x_4(t)$\big) Evolution - Constrained','Interpreter','latex')
xlabel('$x_3(t)$ - [rad]','Interpreter','latex')
ylabel('$x_3(t)$ - [rad/s]','Interpreter','latex')
grid

figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t(1:end-1),us,'-g',LineWidth=1.5)
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]','$u_{max}$','$-u_{max}$','Interpreter', 'latex', 'Location', 'best')
grid

sgtitle('Polytopic Descrption - Constrained Input','Interpreter','latex')


x_store{end+1}=xs;
u_store{end+1}=us;
t_store{end+1}=t;

sgtitle('Polytopic Descrption - Constrained Input','Interpreter','latex')

input('Comparison - Unconstrained Constrained')
close all
clear F t us xs Qe Q1 Q2

% COMPARISON (Unconstrained vs Constrained)
figure(1)
subplot(2,2,1);
plot(t_store{1},x_store{1}(1,:),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(1,:),'g','LineWidth',1.5)
hold off
legend('$x_1(t)$','$x_1(t)$ constrained','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [m]','Interpreter','latex')
grid
subplot(2,2,2);
plot(t_store{1},x_store{1}(2,:),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(2,:),'g','LineWidth',1.5)
hold off
legend('$x_2(t)$','$x_2(t)$ constrained','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Velocity [m/s]','Interpreter','latex')
grid

subplot(2,2,[3 4])
e1.bestPlot([0.5 0.5 0.5])
hold on
e1c.bestPlot([1 0 0])
plot(x_store{1}(1,:),x_store{1}(2,:),'b','LineWidth',1.5);
plot(x_store{2}(1,:),x_store{2}(2,:),'g','LineWidth',1.5);
plot(param.x_eq(1),param.x_eq(2),'ok','LineWidth',1)
hold off
xlabel('$x_1(t)$','Interpreter','latex')
ylabel('$x_2(t)$','Interpreter','latex')
grid
legend('Uncontrained','Constrained','State Evolution Unconstrained','State Evolution Constrained','Equilibrium Point','Interpreter','latex','Location','best')
sgtitle('\big($x_1(t)$,$x_2(t)$\big) Evolution - Unconstrained vs Constrained','Interpreter','latex')

figure(2)
subplot(2,2,1);
plot(t_store{1},x_store{1}(3,:),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(3,:),'g','LineWidth',1.5)
hold off
legend('$x_3(t)$','$x_3(t)$ constrained','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [rad]','Interpreter','latex')
grid
subplot(2,2,2);
plot(t_store{1},x_store{1}(4,:),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(4,:),'g','LineWidth',1.5)
hold off
legend('$x_4(t)$','$x_4(t)$ constrained','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Angular Velocity [rad/s]','Interpreter','latex')
grid

subplot(2,2,[3 4])
e2.bestPlot([0.5 0.5 0.5])
hold on
e2c.bestPlot([1 0 0])
plot(x_store{1}(3,:),x_store{1}(4,:),'b','LineWidth',1.5);
plot(x_store{2}(3,:),x_store{2}(4,:),'g','LineWidth',1.5);
plot(param.x_eq(3),param.x_eq(4),'ok','LineWidth',1)
xlabel('$x_3(t)$','Interpreter','latex')
ylabel('$x_4(t)$','Interpreter','latex')
grid
legend('Uncontrained','Constrained','State Evolution Unconstrained','State Evolution Constrained','Equilibrium Point','Interpreter','latex')
sgtitle('\big($x_3(t)$,$x_4(t)$\big) Evolution - Unconstrained vs Constrained','Interpreter','latex')

figure(3)
plot(t_store{1}(1:end-1),u_store{1},'b','LineWidth',1.5);
hold on
plot(t_store{2}(1:end-1),u_store{2},'g','LineWidth',1.5);
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input - $u(t)$', 'Interpreter', 'latex')
ylim([-5 5])
legend('$u(t)$','$u(t)$ constrained','$u_{max}$','$-u_{max}$','Interpreter','latex','Location','best')
grid
sgtitle('Input Unconstrained vs Input Constrained','Interpreter','latex')

input('Component Wise Constraints')
close all

% POLYTOPIC DESCRIPTION (Component Wise Constraints)
y_max = [0.4, 0.3, 0.4, 0.35]';

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

figure;
subplot(2,2,1)
plot(t,xs(1,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(1),'--k')
yline(y_max(1),'--k')
ylim([-0.5 0.5])
legend('$x_1(t)$ evolution','$-y_{max,1}$','$y_{max,1}$','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [m]','Interpreter','latex')
grid

subplot(2,2,2)
plot(t,xs(2,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(2),'--k')
yline(y_max(2),'--k')
ylim([-0.5 0.5])
legend('$x_2(t)$ evolution','$-y_{max,2}$','$y_{max,2}$','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Velocity [m/s]','Interpreter','latex')
grid

subplot(2,2,3)
plot(t,xs(3,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(3),'--k')
yline(y_max(3),'--k')
ylim([-0.5 0.5])
legend('$x_3(t)$ evolution','$-y_{max,3}$','$y_{max,3}$','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Position [rad]','Interpreter','latex')
grid

subplot(2,2,4)
plot(t,xs(4,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(4),'--k')
yline(y_max(4),'--k')
ylim([-0.5 0.5])
legend('$x_4(t)$ evolution','$-y_{max,4}$','$y_{max,4}$','Interpreter','latex')
xlabel('Time [s]','Interpreter','latex')
ylabel('Angular Velocity [rad/s]','Interpreter','latex')
grid

sgtitle('State Evolution under Constrains','Interpreter','latex')

figure
plot(t(1:end-1),us,'-g',LineWidth=1.5)
hold on
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]','Interpreter','latex')
ylabel('Input - $u(t)$','Interpreter','latex')
ylim([-4 4])
legend('$u(t)$ - Applied Torque [rad/s$^2$]','$u_{max}$','-$u_{max}$','Interpreter','latex')
grid

sgtitle('Input $u(t)$ Constrained under State Constrains','Interpreter','latex')

% ELLIPSOID (Component Wise Constraints)

Q1 = QQ(1:2,1:2);
Q2 = QQ(3:4,3:4);

e1cc = ell(Q1);
e2cc = ell(Q2);

figure
subplot(2,1,1)
e1.bestPlot([0.7 0.7 0.7],1)
hold on
e1cc.bestPlot([1 0 0])
plot(x_store{1}(1,:),x_store{1}(2,:),'b','LineWidth',1.5);
plot(xs(1,:),xs(2,:),'g','LineWidth',1.5);
plot(param.x_eq(1),param.x_eq(2),'ok','LineWidth',1)
rectangle('Position',[-y_max(1),-y_max(2), 2*y_max(1),2*y_max(2)],'EdgeColor', [1 0.5 0],'LineStyle','--','LineWidth',1.5)
hold off
xlabel('$x_1(t)$','Interpreter','latex')
ylabel('$x_2(t)$','Interpreter','latex')
grid
legend('Uncontrained','Constrained','State Evolution Unconstrained','State Evolution Constrained','Equilibrium Point','Interpreter','latex')
title('\big($x_1(t)$,$x_2(t)$\big) Evolution - Constrained Input and State','Interpreter','latex')

subplot(2,1,2)
e2.bestPlot([0.7 0.7 0.7],1)
hold on
e2cc.bestPlot([1 0 0])
plot(x_store{1}(3,:),x_store{1}(4,:),'b','LineWidth',1.5);
plot(xs(3,:),xs(4,:),'g','LineWidth',1.5);
plot(param.x_eq(1),param.x_eq(2),'ok','LineWidth',1)
rectangle('Position',[-y_max(3),-y_max(4), 2*y_max(3),2*y_max(4)],'EdgeColor', [1 0.5 0],'LineStyle','--','LineWidth',1.5)
hold off
xlabel('$x_3(t)$','Interpreter','latex')
ylabel('$x_4(t)$','Interpreter','latex')
grid
legend('Uncontrained','Constrained','State Evolution Unconstrained','State Evolution Constrained','Equilibrium Point','Interpreter','latex')
title('\big($x_3(t)$,$x_4(t)$\big) Evolution - Constrained Input and State','Interpreter','latex')

