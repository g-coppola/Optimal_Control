classdef MPC
    properties
        A, B, C
        f
        ff
        param
    end
    
    methods (Access = public)
        function obj = MPC(f,ft,param,A,B,C)
            obj.f = f;
            obj.ff = ft;
            obj.param = param;
            obj.A = A;
            obj.B = B;
            obj.C = C;
        end
        
        function [x, u, t, Qe] = unconstrainedRHC(obj,x0)
            [F,Qe] = obj.uncRHC(x0);
            
            [t,x] = ode45(@(t,x)obj.f(x-obj.param.x_eq, obj.param.u_eq + F*(x-obj.param.x_eq)), obj.param.tspan, x0);
            
            u = zeros(size(t));
            for i = 1:length(t)
                u(i) = obj.param.u_eq + F * (x(i,:)' - obj.param.x_eq);
            end
        end

        function [x, u, t, Qe] = constrainedInputRHC(obj,x0,u_max)
            [F,Qe] = obj.inConRHC(x0,u_max);
            
            [t,x] = ode45(@(t,x)obj.f(x-obj.param.x_eq, obj.param.u_eq + F*(x-obj.param.x_eq)), obj.param.tspan, x0);
            
            u = zeros(size(t));
            for i = 1:length(t)
                u(i) = obj.param.u_eq + F * (x(i,:)' - obj.param.x_eq);
            end
        end

        function [x, u, t, Qe] = compWiseOutputRHC(obj,x0,u_max,y_max)
            [F,Qe] = obj.outCompWiseRHC(x0,u_max,y_max);
            
            [t,x] = ode45(@(t,x)obj.f(x-obj.param.x_eq, obj.param.u_eq + F*(x-obj.param.x_eq)), obj.param.tspan, x0);
            
            u = zeros(size(t));
            for i = 1:length(t)
                u(i) = obj.param.u_eq + F * (x(i,:)' - obj.param.x_eq);
            end
        end

        function [F, P] = offLine_phase(obj,xN,u_max,y_max)
            [F,Qe] = obj.outCompWiseRHC(xN,u_max,y_max);
            P = inv(Qe);
        end

        function [x, u, t, tu] = onLine_phase(obj,ft,x0,u_max,y_max,Q,R,N,h,F,P)
            n = obj.param.n;
            m = obj.param.m;
            x = zeros(n,N+1);
            x(:,1) = x0;
            xk = x0;

            lb_X = repmat(-y_max, h + 1, 1);
            ub_X = repmat(y_max, h + 1, 1);
            
            lb_U = repmat(-u_max * ones(m, 1), h, 1);
            ub_U = repmat(u_max * ones(m, 1), h, 1);
            
            lb = [lb_X; lb_U]; 
            ub = [ub_X; ub_U];

            Aineq = [];
            bineq = [];
            Aeq = [];
            beq = [];
            
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
            
            U0 = zeros(m,h);
            X0 = zeros(n,h+1);
            X0(:,1) = xk;

            x_sim = xk;

            for k = 1:h
                u_sim = U0(:, k);
                dx = ft(0, x_sim, u_sim);
                x_sim = x_sim + obj.param.dt*dx;
                X0(:, k+1) = x_sim;
            end

            Z0 = [X0(:); U0(:)];

            u = zeros(m,N);

            for k = 1:N
                if k>h
                    uk = F*xk; 
                else
                    [Z_ottimo, ~, exitflag, ~] = fmincon(@(Z)obj.objFunct(Z,h,Q,R,P),Z0,Aineq, bineq, Aeq, beq, lb, ub, @(Z)obj.nonLinConstraints(Z,xk,h,P,ft),options);

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
                         uk = zeros(m,1);
                    end
                end

                u(:,k) = uk;
                dx = ft(0,xk,uk);
                x(:,k+1) = xk+obj.param.dt*dx;
                
                xk = x(:,k+1);
                Z0(1:n) = xk;
            end

            t= 0:obj.param.dt:(N*obj.param.dt);
            tu = 0:obj.param.dt:((N-1)*obj.param.dt);

        end

    end

    methods (Access = private)
        function [F,Qe] = uncRHC(obj,x0)
            n = obj.param.n;
            m = obj.param.m;
            Ac = obj.A;
            Bc = obj.B;

            yalmip('clear')
            gamma=sdpvar;
            Q=sdpvar(n);
            Y=sdpvar(m,n,'full');

            R = obj.param.RR;
            Q1 = obj.param.Q1;

            con=[Q>=0, [1 x0';x0 Q]>=0];
            for j=1:length(Ac)
                Aj=Ac{j};
                Bj=Bc{j};
                con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                            Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                            sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                            sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0];
            end

            % Define solver options
            options = sdpsettings('solver', 'sedumi', 'verbose', 0);

            % Solve the optimization problem
            sol = optimize(con, gamma, options);
            
            % Check if the problem is solved
            if sol.problem == 0
                Qe = inv(value(Q));
                F=value(Y)*Qe;
            else
                disp(sol.info);
                error('Something went wrong:');
            end
        end

        function [F, Qe] = inConRHC(obj,x0,u_max)
            n = obj.param.n;
            m = obj.param.m;
            Ac = obj.A;
            Bc = obj.B;

            R = obj.param.RR;
            Q1 = obj.param.Q1;

            yalmip('clear')
            gamma=sdpvar;
            Q=sdpvar(n);
            Y=sdpvar(m,n,'full');

            con=[Q>=1e-6*eye(n), [1 x0';x0 Q]>=0];
            for j=1:length(Ac)
                Aj=Ac{j};
                Bj=Bc{j};
                con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                            Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                            sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                            sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0];
            end
            
            con = [con, [u_max^2*eye(m) Y; Y' Q]>=0];
            
            % Define solver options
            options = sdpsettings('solver', 'sedumi', 'verbose', 0, 'debug', 1);
            % Solve the optimization problem
            sol = optimize(con, gamma, options);
            % Check if the problem is solved
            if sol.problem == 0
                Qe = inv(value(Q));
                F=value(Y)*Qe;
            else
                disp(sol.info);
                error('Something went wrong:');
            end
        end

        function [F,Qe] = outCompWiseRHC(obj,x0,u_max,y_max)
            n = obj.param.n;
            m = obj.param.m;
            ny = size(obj.C,1);
            
            Ac = obj.A;
            Bc = obj.B;

            R = obj.param.RR;
            Q1 = obj.param.Q1;

            yalmip('clear')
            Q=sdpvar(n);
            gamma=sdpvar;
            Y=sdpvar(m,n,'full');

            con=[Q>=0, [1 x0';x0 Q]>=0];
            for j = 1:length(Ac)
                Aj = Ac{j};
                Bj = Bc{j};
                con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                            Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                            sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                            sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0];
                for l = 1:ny
                    Cl = obj.C(l,:);
                    yl = y_max(l);
                    nc = size(Cl,1);
                    con = [con, [Q (Aj*Q+Bj*Y)'*Cl'; Cl*(Aj*Q+Bj*Y) (yl^2)*eye(nc)]>=0];
                end
            end
            
            con = [con, [u_max^2*eye(m) Y; Y' Q]>=0];
            
            % Define solver options
            options = sdpsettings('solver', 'sedumi', 'verbose', 0);
            % Solve the optimization problem
            sol = optimize(con, gamma, options);
            % Check if the problem is solved
            if sol.problem == 0
                Qe = inv(value(Q));
                F=value(Y)*Qe;
            else
                disp(sol.info);
                error('Something went wrong:');
            end 
        end

        function [J] = objFunct(obj,Z,h,Q,R,P)
            J = 0;
            Ac = obj.A;
            Bc = obj.B;
            n = obj.param.n;
            m = obj.param.m;
            l = size(Ac,2);
            
            num_X_vars = n * (h + 1);
            X_vec = Z(1:num_X_vars);
            U_vec = Z(num_X_vars + 1 : end);
            
            X = reshape(X_vec, n, h + 1);
            U = reshape(U_vec, m, h);
         
            for k = 1:h
                u = U(:,k); 
                x = X(:,k);
                c = zeros(1,l);
                for i = 1:l
                    c(i) = (Ac{i}*x+Bc{i}*u)'*Q*(Ac{i}*x+Bc{i}*u)+ u'*R *u;
                end
                J = J + max(c);
            end
        
            % Terminal Cost
            cf = zeros(1,l);
            for i = 1:l
                cf(i) = (Ac{i}*X(:,h)+Bc{i}*U(:,h))'*P*(Ac{i}*X(:,h)+Bc{i}*U(:,h));
            end
            J = J + max(cf);
        end
    
        function [c,ceq] = nonLinConstraints(obj,Z,xk,h,P,ft)
            c = [];
            ceq = [];
            dt = obj.param.dt;
            Ac = obj.A;
            Bc = obj.B;
        
            n = obj.param.n;
            m = obj.param.m;
            l = size(Ac,2);
            
            num_X_vars = n * (h + 1);
            X_vec = Z(1:num_X_vars);
            U_vec = Z(num_X_vars + 1 : end);
            
            X = reshape(X_vec, n, h + 1);
            U = reshape(U_vec, m, h);
            
            for i = 1:l
                cf = (Ac{i}*X(:,h)+Bc{i}+U(:,h))'*P*(Ac{i}*X(:,h)+Bc{i}+U(:,h))-1;
                c = [c; cf];
            end
        
            X0 = X(:,1);
            c0 = X0-xk;
            ceq = [ceq; c0]; 
        
            for k = 1:h
                X_pre = X(:,k);
                X_post = X(:,k+1);
                Uk = U(:,k); 
                
                dx = ft(0,X_pre,Uk);
                Xk_post = X_pre+dt*dx; 
        
                cc = X_post - Xk_post;
                ceq = [ceq; cc];
            end
end
    end

end

