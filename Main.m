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
param.r = 0.12;

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
    -(m*g)*((x(1)*cos(x(3)))/(I+m*x(1)^2))-2*m*((x(1)*x(2)*x(4))/(I+m*x(1)^2))+u/(I+m*x(1)^2)
    ];

% Equilibrium
param.x_eq = [0 0 0 0]';
param.u_eq = 0;

% Initial Conditions
x0 = param.x_eq + [0.1 0 -0.05 0]';

%% FREE RESPONSE
f_fun = matlabFunction(f,'Vars',{x,u});

[t,xs] = ode45(@(t,x)f_fun(x,0),param.tspan,x0);

figure;
plot(t,xs,LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
title('Free Response','Interpreter','latex')
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
grid

input('Press to LQR (Zero Regulation)')
close all   
clear xs t

%% LQR (Zero Regulation)
N = round(param.tspan(2)/param.dt);

ft = @(t, x, u) [x(2);
    -(m*g/M)*sin(x(3))+(m/M)*x(1)*x(4)^2;
    x(4);
    -(m*g)*((x(1)*cos(x(3)))/(I+m*x(1)^2))-2*m*((x(1)*x(2)*x(4))/(I+m*x(1)^2))+u/(I+m*x(1)^2)
    ];

param.Q = diag([10 1 5 1]);
param.R = 0.01*eye(param.m);
param.S = eye(param.n);

OC = OptControl(param,f,x,u);

% Control for infinite Horizon
[t, xs, us, e] = OC.lqr(ft,x0,N);

figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('States $\mathbf{x}(t)$', 'Interpreter', 'latex')
title('LQR (Infinite) - Zero Regulation', 'Interpreter', 'latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t,us,'-g',LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]', 'Interpreter', 'latex', 'Location', 'best')
grid

input("")
clear xs us t

% Predictive Control for Finite Horizon
[t, xs,us] = OC.OLQR(ft,x0,N);

figure;
subplot(2,1,1)
plot(t,xs,LineWidth=1.5);
xlabel('Time [s]')
ylabel('States - x(t)')
title('LQR (Model Predictive Finite) - Zero Regulation', 'Interpreter','latex')
grid
legend('$x_1(t)$ - Ball Position [m]', '$x_2(t)$ - Ball Velocity [m/s]', '$x_3(t)$ - Beam Angle [rad]', '$x_4(t)$ - Beam Velocity [rad/s]', 'Interpreter', 'latex', 'Location', 'best')
hold off
subplot(2,1,2)
plot(t,us,'-g',LineWidth=1.5)
xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('Input $u(t)$', 'Interpreter', 'latex')
legend('$u(t)$ - Applied Torque [rad/s$^2$]', 'Interpreter', 'latex', 'Location', 'best')
grid

input('Press to Polytopic Description (Unconstrained)')
close all
clear xs us t

%% POLYTOPIC DESCRIPTION (Unconstrained)
L = 0.5;
th = pi/8;

tspan = [0 15];

x_store = {};
t_store = {};

[A_cell,B_cell] = build_pol(param,dt,L,th);

[F,Qe] = uncostrainedRHC(A_cell,B_cell,x0);

[t,x] = ode45(@(t,x)f_fun(x-x_eq, u_eq + F*(x-x_eq)), tspan, x0);

u = zeros(size(t));
for i = 1:length(t)
    u(i) = u_eq + F * (x(i,:)' - x_eq);
end

x_store{end+1}=x;
t_store{end+1}=t;

figure;
subplot(2,1,1)
plot(t,x,LineWidth=1.5);
xlabel('Time [s]')
ylabel('States - x(t)')
title('Polytopic Descrption - Unconstrained')
grid
legend('x_1(t) - Ball Position [m]', 'x_2(t) - Ball Velocity [m/s]', 'x_3(t) - Beam Angle [rad]','x_4(t) - Beam Velocity [rad/s]')
hold off
subplot(2,1,2)
plot(t,u,'-g',LineWidth=1.5)
xlabel('Time [s]')
ylabel('Input - u(t)')
legend('u(t) - Applyed Torque [rad/s^2]')
grid

input('')

%% ELLYPSOID (unconstrained)
Q1 = Qe(1:2,1:2);
Q2 = Qe(3:4,3:4);

e1_unc = ell(Q1);
e2_unc = ell(Q2);

figure(2)
subplot(2,1,1)
e1_unc.bestPlot([1 0 0])
hold on
plot(x(:,1),x(:,2),'b','LineWidth',1.5);
xlabel('x_1 - [m]')
ylabel('x_2 - [m/s]')
legend('Ellipsoid','Equilibrium Point','State Evolution')
title('(x_1,x_2) Evolution - Uncontrained')

subplot(2,1,2)
e2_unc.bestPlot([1 0 0])
hold on
plot(x(:,3),x(:,4),'b','LineWidth',1.5);
xlabel('x_3 - [rad]')
ylabel('x_4 - [rad/s]')
legend('Ellipsoid','Equilibrium Point','State Evolution')
title('(x_3,x_4) Evolution - Uncontrained')

input('Press to Polytopic Description (Constrained Input)')
close all
clear F t u x Qe Q1 Q2

% POLYTOPIC DESCRIPTION (Constrained Input)
u_max = 4.5 + u_eq;

[F,Qe] = constrainedRHC(A_cell,B_cell,x0,u_max);

[t,x] = ode45(@(t,x)f_fun(x-x_eq,u_eq+F*(x-x_eq)),tspan,x0);

u = zeros(size(t));
for i = 1:length(t)
    u(i) = u_eq + F * (x(i,:)' - x_eq);
end

x_store{end+1}=x;
t_store{end+1}=t;

figure;
subplot(2,1,1)
plot(t,x,LineWidth=1.5);
xlabel('Time [s]')
ylabel('States - x(t)')
title('Polytopic Descrption - Constrained Input')
grid
legend('x_1(t) - Ball Position [m]', 'x_2(t) - Ball Velocity [m/s]', 'x_3(t) - Beam Angle [rad]','x_4(t) - Beam Velocity [rad/s]')
hold off
subplot(2,1,2)
plot(t,u,'-g',LineWidth=1.5)
hold on
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]')
ylabel('Input - u(t)')
ylim([-7.2 7.2])
legend('u(t) - Applyed Torque [rad/s^2]','u_{max}','+u_{max}')
grid

input('')

% ELYPSOID (Constrained Input)
Q1 = Qe(1:2,1:2);
Q2 = Qe(3:4,3:4);

e1_con = ell(Q1);
e2_con = ell(Q2);

figure(2)
subplot(2,1,1)
e1_con.bestPlot([1 0 0])
hold on
plot(x(:,1),x(:,2),'b','LineWidth',1.5);
xlabel('x_1 - [m]')
ylabel('x_2 - [m/s]')
legend('Ellipsoid','Equilibrium Point','State Evolution')
title('(x_1,x_2) Evolution - Contrained Input')

subplot(2,1,2)
e2_con.bestPlot([1 0 0])
hold on
plot(x(:,3),x(:,4),'b','LineWidth',1.5);
xlabel('x_3 - [rad]')
ylabel('x_4 - [rad/s]')
legend('Ellipsoid','Equilibrium Point','State Evolution')
title('(x_3,x_4) Evolution - Contrained Input')

input('Comparison - Unconstrained Constrained')
close all
clear F t u x Qe Q1 Q2

% COMPARISON (Unconstrained vs Constrained)
figure(1)
subplot(2,2,1);
plot(t_store{1},x_store{1}(:,1),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(:,1),'g','LineWidth',1.5)
hold off
legend('x_1','x_1 constrained')
xlabel('Time [s]')
ylabel('Position [m]')
grid
subplot(2,2,2);
plot(t_store{1},x_store{1}(:,2),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(:,2),'g','LineWidth',1.5)
hold off
legend('x_2','x_2 constrained')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid

subplot(2,2,[3 4])
e1_unc.bestPlot([0.5 0.5 0.5])
hold on
e1_con.bestPlot([1 0 0])
plot(x_store{1}(:,1),x_store{1}(:,2),'b','LineWidth',1.5);
plot(x_store{2}(:,1),x_store{2}(:,2),'g','LineWidth',1.5);
xlabel('x_1')
ylabel('x_2')
grid
legend('Uncontrained','Equilibrium Point','Constrained','Equilibrium Point','State Evolution Unconstrained','State Evolution Constrained')
sgtitle('(x_1,x_2) evolution - Unconstrained vs Constrained Input')

figure(2)
subplot(2,2,1);
plot(t_store{1},x_store{1}(:,3),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(:,3),'g','LineWidth',1.5)
hold off
legend('x_3','x_3 constrained')
xlabel('Time [s]')
ylabel('Position [rad]')
grid
subplot(2,2,2);
plot(t_store{1},x_store{1}(:,4),'b','LineWidth',1.5)
hold on
plot(t_store{2},x_store{2}(:,4),'g','LineWidth',1.5)
hold off
legend('x_4','x_4 constrained')
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
grid

subplot(2,2,[3 4])
e2_unc.bestPlot([0.5 0.5 0.5])
hold on
e2_con.bestPlot([1 0 0])
plot(x_store{1}(:,3),x_store{1}(:,4),'b','LineWidth',1.5);
plot(x_store{2}(:,3),x_store{2}(:,4),'g','LineWidth',1.5);
xlabel('x_3')
ylabel('x_4')
grid
legend('Uncontrained','Equilibrium Point','Constrained','Equilibrium Point','State Evolution Unconstrained','State Evolution Constrained')
sgtitle('(x_3,x_4) evolution - Unconstrained vs Constrained Input')

input('Component Wise Constraints')
close all

% POLYTOPIC DESCRIPTION (Component Wise Constraints)
y_max = [0.2 0.32 0.18 0.76]';
[F,Qe] = compWiseCon(A_cell,B_cell,C,x0,y_max,u_max);

[t,x] = ode45(@(t,x)f_fun(x-x_eq,u_eq+F*(x-x_eq)),tspan,x0);

u = zeros(size(t));
for i = 1:length(t)
    u(i) = u_eq + F * (x(i,:)' - x_eq);
end

figure;
subplot(2,2,1)
plot(t,x(:,1),'-b',LineWidth=1.5);
hold on
yline(-y_max(1),'--k')
yline(y_max(1),'--k')
ylim([-0.32 0.32])
legend('x_1 evolution','- y_{max,1}','y_{max,1}')
xlabel('Time [s]')
ylabel('Position [m]')
grid

subplot(2,2,2)
plot(t,x(:,2),'-b',LineWidth=1.5);
hold on
yline(-y_max(2),'--k')
yline(y_max(2),'--k')
ylim([-0.7 0.7])
legend('x_2 evolution','- y_{max,2}','y_{max,2}')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid

subplot(2,2,3)
plot(t,x(:,3),'-b',LineWidth=1.5);
hold on
yline(-y_max(3),'--k')
yline(y_max(3),'--k')
ylim([-0.31 0.31])
legend('x_3 evolution','- y_{max,3}','y_{max,3}')
xlabel('Time [s]')
ylabel('Angle [rad]')
grid

subplot(2,2,4)
plot(t,x(:,4),'-b',LineWidth=1.5);
hold on
yline(-y_max(4),'--k')
yline(y_max(4),'--k')
ylim([-1 1])
legend('x_4 evolution','- y_{max,4}','y_{max,4}')
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
grid

% ELLIPSOID (Component Wise Constraints)

Q1 = Qe(1:2,1:2);
Q2 = Qe(3:4,3:4);

e1_ccon = ell(Q1);
e2_ccon = ell(Q2);

figure
plot(t,u,'-g',LineWidth=1.5)
hold on
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]')
ylabel('Input - u(t)')
ylim([-7.2 7.2])
legend('u(t) - Applyed Torque [rad/s^2]','u_{max}','+u_{max}')
grid

figure
subplot(2,1,1)
e1_unc.bestPlot([0.7 0.7 0.7],1)
hold on
e1_ccon.bestPlot([1 0 0])
plot(x_store{1}(:,1),x_store{1}(:,2),'b','LineWidth',1.5);
plot(x(:,1),x(:,2),'g','LineWidth',1.5);
rectangle('Position',[-y_max(1),-y_max(2), 2*y_max(1),2*y_max(2)],'EdgeColor', [1 0.5 0],'LineStyle','--','LineWidth',1.5)
xlabel('x_1')
ylabel('x_2')
grid
legend('Uncontrained','Equilibrium Point','Constrained','Equilibrium Point','State Evolution Unconstrained','State Evolution Constrained')

subplot(2,1,2)
e2_unc.bestPlot([0.7 0.7 0.7],1)
hold on
e2_ccon.bestPlot([1 0 0])
plot(x_store{1}(:,3),x_store{1}(:,4),'b','LineWidth',1.5);
plot(x(:,3),x(:,4),'g','LineWidth',1.5);
rectangle('Position',[-y_max(3),-y_max(4), 2*y_max(3),2*y_max(4)],'EdgeColor', [1 0.5 0],'LineStyle','--','LineWidth',1.5)

xlabel('x_3')
ylabel('x_4')
grid
legend('Uncontrained','Equilibrium Point','Constrained','Equilibrium Point','State Evolution Unconstrained','State Evolution Constrained')

input('Online Phase!!')
close all;

%% ONLINE MPC SIMULATION - DUAL MODE CONTROLLER

N = 80;
h = 12;

P = inv(Qe);

Q = diag([10 10 0.1 0.1]);
R = 0.1;

u_max = 10;
y_max = [0.8 1.5 pi/4 1.7]'; 

x_online = zeros(param.n,N+1);
x0 = [0.6, 0, pi/8, 0]';
x_online(:,1) = x0;
xk = x0; 


lb_X = repmat(-y_max, h + 1, 1);
ub_X = repmat(y_max, h + 1, 1);

lb_U = repmat(-u_max * ones(param.m, 1), h, 1);
ub_U = repmat(u_max * ones(param.m, 1), h, 1);

lb = [lb_X; lb_U]; 
ub = [ub_X; ub_U];

Aineq = [];
bineq = [];
Aeq = [];
beq = [];

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

U0 = zeros(param.m,h);
X0 = zeros(param.n,h+1);
X0(:,1) = xk;

x_sim = xk;
for k = 1:h
    u_sim = U0(:, k);
    dx = ft(0, x_sim, u_sim);
    x_sim = x_sim + dt*dx;
    X0(:, k+1) = x_sim;
end

Z0 = [X0(:); U0(:)];

u_online = zeros(param.m, N);

n = param.n; 
m = param.m;

%ccc = zeros(2,h);

for k = 1:N
    %ccc(:,k) = [e1_ccon.isInternal([xk(1);xk(2)]); e2_ccon.isInternal([xk(3);xk(4)])]
   
    if k>h
        uk = F*xk;
    else
          
    [Z_ottimo, ~, exitflag, ~] = fmincon(@(Z)objFunct(Z,h,A_cell,B_cell,Q,R,P,param),Z0,Aineq, bineq, Aeq, beq, lb, ub,@(Z)nonLinConstraints(Z,xk,A_cell,B_cell,h,P,dt,ft,param),options);
    
    uk = zeros(m,1); 
    
    if exitflag > 0
        num_X_vars = n * (h + 1);
        
        
        X_opt_vec = Z_ottimo(1:num_X_vars);
        U_opt_vec = Z_ottimo(num_X_vars + 1 : end);
        
        X_opt = reshape(X_opt_vec, n, h + 1);
        U_opt = reshape(U_opt_vec, m, h);

        uk = U_opt(:, 1); 

        x_opt_pr = [X_opt(:, 2:end), X_opt(:, end)];
        u_opt_pr = [U_opt(:, 2:end), U_opt(:, end)];

        Z0 = [x_opt_pr(:); u_opt_pr(:)];
        
    else
        fprintf('Attenzione: fmincon non ha trovato una soluzione ottimale alla iterazione %d\n', k);
        uk = zeros(m,1);
    end
    end
    
    u_online(:,k) = uk;
    
    dx = ft(0,xk,uk);
    x_online(:,k+1) = xk+dt*dx;

    xk = x_online(:,k+1); 

    Z0(1:n) = xk;
end


figure
subplot(2,1,1)
hold on
e1_ccon.bestPlot([1 0 0]);
plot(x_online(1,:),x_online(2,:),'g','LineWidth',1.5);
rectangle('Position',[-y_max(1),-y_max(2), 2*y_max(1),2*y_max(2)],'EdgeColor', [1 0.5 0],'LineStyle','--','LineWidth',1.5)
xlabel('x_1')
ylabel('x_2')
legend('Ellipsoid','Equilibrium Point','State Evolution')
grid on

subplot(2,1,2)
hold on
e2_ccon.bestPlot([1 0 0]);
plot(x_online(3,:),x_online(4,:),'g','LineWidth',1.5);
rectangle('Position', [-y_max(3), -y_max(4), 2*y_max(3), 2*y_max(4)],'EdgeColor', [1 0.5 0],'LineStyle', '--', 'LineWidth', 1.5)
xlabel('x_3')
ylabel('x_4')
legend('Ellipsoid','Equilibrium Point','State Evolution')
grid on

t_online = 0:dt:(N*dt);
tu = 0:dt:((N-1)*dt);

figure;
subplot(2,2,1)
plot(t_online,x_online(1,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(1),'--k')
yline(y_max(1),'--k')
ylim([-1.7 1.7])
legend('x_1 evolution','- y_{max,1}','y_{max,1}')
xlabel('Time [s]')
ylabel('Position [m]')
grid

subplot(2,2,2)
plot(t_online,x_online(2,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(2),'--k')
yline(y_max(2),'--k')
ylim([-1.9 1.9])
legend('x_2 evolution','- y_{max,2}','y_{max,2}')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid

subplot(2,2,3)
plot(t_online,x_online(3,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(3),'--k')
yline(y_max(3),'--k')
ylim([-1.5 1.5])
legend('x_3 evolution','- y_{max,3}','y_{max,3}')
xlabel('Time [s]')
ylabel('Angle [rad]')
grid

subplot(2,2,4)
plot(t_online,x_online(4,:),'-b',LineWidth=1.5);
hold on
yline(-y_max(4),'--k')
yline(y_max(4),'--k')
ylim([-1.9 1.9])
legend('x_4 evolution','- y_{max,4}','y_{max,4}')
xlabel('Time [s]')
ylabel('Angular Velocity [rad/s]')
grid

sgtitle('x(t) - States Evolutions')

figure
plot(tu,u_online,'-g',LineWidth=1.5)
hold on
yline(u_max,'--k');
yline(-u_max,'--k')
xlabel('Time [s]')
ylabel('Input - u(t)')
ylim([-11.5 11.5])
legend('u(t)','u_{max}','+u_{max}')
title('u(t) - Applyed Torque [rad/s^2]')
grid
