classdef ell
    
    properties
        Q
        q
        n
    end
    
    methods
        function obj = ell(Q,q)
            obj.Q = Q;
            obj.n = size(Q,1);
            
            if nargin == 1 % se il numero di input è 1 allora q = 0n
                q = zeros(obj.n,1);
            end
            obj.q = q;
        end
        
        
        function plot(obj)
            yalmip("clear");
            x = sdpvar(2,1);
            con = [1 (x-obj.q)'; (x-obj.q) obj.Q] >= 0;
            figure
            plot(con);
            hold on
            plot(obj.q(1),obj.q(2),'ob','LineWidth',1.5)
            hold off
            grid
        end

        function b = isInternal(obj,x)
            b = (x-obj.q)'*obj.Q*(x-obj.q) <= 1;
        end

        function [phi,b] = supp(obj,z)
            phi = z'*obj.q+sqrt(z'*inv(obj.Q)*z);
            b = inv(obj.Q)*z/sqrt(z'*inv(obj.Q)*z)+obj.q;
        end

        function bestPlot(obj,color,bool)
            if nargin == 2 
                bool = 0;
            end
            % Plot con Support Function
            delta = 0.01;
            Theta = 0:delta:2*pi+delta;
            
            x = zeros(size(obj.q,1),length(Theta));
            
            for i = 1:length(Theta)
                z = [cos(Theta(i)) sin(Theta(i))]';
                [~, b] = obj.supp(z);
                x(:,i) = b;
            end

            if bool == 0
                plot(x(1,:),x(2,:),'Color',color,'LineWidth',1)
                hold on
                plot(obj.q(1),obj.q(2),'ob','LineWidth',1)
                grid
            else
                plot(x(1,:),x(2,:),'Color',color,'LineWidth',1,'LineStyle', '--')
                hold on
                plot(obj.q(1),obj.q(2),'ob','LineWidth',1)
                grid
            end
            
        end
    end
end

