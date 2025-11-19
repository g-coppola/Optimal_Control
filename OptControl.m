classdef OptControl
    properties
        A, B
        Q, R, S, SS, PP
        param
        f, xsym, usym
    end
    
    methods (Access = public)
        function obj = OptControl(param,f,state,input)           
            obj.Q = param.Q;
            obj.R = param.R;
            obj.S = param.S;
            obj.f = f;
            obj.param = param;
            obj.xsym = state;
            obj.usym = input;
            [obj.A, obj.B] = obj.linearize();
        end
        
        function [t, x, u, e] = lqr(obj, ft, x0,N)
            [K, obj.SS, e] = dlqr(obj.A,obj.B,obj.Q,obj.R);

            x = zeros(obj.param.n,N);
            u = zeros(obj.param.m,N);
            x(:,1) = x0;

            for k = 1:N-1
                u(:,k) = obj.param.u_eq-K*(x(:,k)-obj.param.x_eq);
                dx = ft(0,x(:,k),u(:,k));
                x(:,k+1) = x(:,k)+obj.param.dt*dx;
            end
            u(N) = obj.param.u_eq-K*(x(:,N)-obj.param.x_eq);

            t = 0:obj.param.dt:obj.param.tspan(2)-obj.param.dt;
        end

        function [t, x, u] = OLQR(obj,ft,x0,N)
            P = cell(1,N);
            P{N} = obj.S;
            F = cell(1,N);
            
                for k = N-1:-1:1
                    F{k+1} = inv(obj.R + obj.B'*P{k+1}*obj.B)*obj.B'*P{k+1}*obj.A;
                    P{k} = obj.A'*P{k+1}*obj.A + obj.Q - obj.A'*P{k+1}*obj.B*inv(obj.R+obj.B'*P{k+1}*obj.B)*obj.B'*P{k+1}*obj.A;
                end

            F{1} = inv(obj.R + obj.B'*P{2}*obj.B)*obj.B'*P{2}*obj.A;

            obj.PP = P;

            x = zeros(obj.param.n,N);
            u = zeros(obj.param.m,N);
            x(:,1) = x0;

            for k = 1:N-1
                u(:,k) = obj.param.u_eq-F{k}*(x(:,k)-obj.param.x_eq);
                dx = ft(0,x(:,k),u(:,k));
                x(:,k+1) = x(:,k)+obj.param.dt*dx;
            end
            u(N) = obj.param.u_eq-F{N}*(x(:,N)-obj.param.x_eq);

            t = 0:obj.param.dt:obj.param.tspan(2)-obj.param.dt;
        end
    end

    methods (Access = private)
        function [A, B] = linearize(obj)
            A_tilde = jacobian(obj.f,obj.xsym);
            B_tilde = jacobian(obj.f,obj.usym);

            xeq = obj.param.x_eq;
            ueq = obj.param.u_eq;

            A_eq = matlabFunction(A_tilde,'Vars',{obj.xsym,obj.usym});
            B_eq = matlabFunction(B_tilde,'Vars',{obj.xsym,obj.usym});

            Al  = A_eq(xeq,ueq);
            Bl = B_eq(xeq,ueq);
            n = size(Al,1);

            Cl = eye(n);
            Dl = 0;

            dt = obj.param.dt;

            sysc = ss(Al,Bl,Cl,Dl);
            sysd = c2d(sysc, dt, 'zoh');
            A = sysd.A;
            B = sysd.B;
        end
    end
end

