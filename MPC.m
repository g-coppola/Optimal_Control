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
        
        function [F, Qe] = unconstrainedRHC(obj,xk)
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

            con=[Q>=0, [1 xk';xk Q]>=0];
            for j=1:length(Ac)
                Aj=Ac{j};
                Bj=Bc{j};
                con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                            Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                            sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                            sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0];
            end

            % % Define solver options
            % options = sdpsettings('solver', 'sedumi', 'verbose', 0);

            % Solve the optimization problem
            sol = optimize(con, gamma);
            
            % Check if the problem is solved
            if sol.problem == 0
                Qe = inv(value(Q));
                F=value(Y)*Qe;
            else
                disp(sol.info);
                error('Something went wrong:');
            end
        end

        function [F, Qe] = constrainedInputRHC(obj,xk,u_max)
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

            con=[Q>=1e-6*eye(n), [1 xk';xk Q]>=0];
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

        function [F, Qe] = compWiseOutputRHC(obj,xk,u_max,y_max)
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

            con=[Q>=0, [1 xk';xk Q]>=0];
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
    end
end


