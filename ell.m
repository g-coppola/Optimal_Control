classdef ell
    
    properties
        P % Symmetric positive definite matrix
        q % Center of the ellipsoid
        n % Dimension of the ellipsoid
        tol = 1e-07; % Tollerance for the projection
    end
    
    methods
        function obj = ell(P,q)
            %ELL construct an instance of this class
            
            if nargin == 0
                return;
            end

            obj.P = P;
            obj.n = size(P,1);

            if nargin == 1
                q = zeros(obj.n,1);
            end
            obj.q = q;

        end


        function b = isinternal(obj,x,B,flag)

            % This method verifies that the x is inside the ellipsoid
            
            if nargin > 2

                if length(x) > size(B,1) && length(x) == size(B,2)
                    x = B*x;
                end

                checkBvalidity(obj,B,length(x));

                if nargin < 4
                    flag = 2;
                end

                switch(flag)
                    case 1
                        e_to_use = obj.sectEll(B);
                    otherwise
                        e_to_use = obj.projEll(B);
                end

                b = e_to_use.isinternal(x);

            elseif nargin == 2

                b = (x-obj.q)'*obj.P*(x-obj.q) <= 1 + 1e-6; % To avoid numerical problems
            else

                error('x not provided');
            end
        end
        
        function [x, center] = get2DBoundary(obj,B,flag)
            % This method returns the boundary and the center of the
            % ellipsoid

            if obj.n == 2
               
                x = [];
                res = 150;
                theta = 0:pi/res:2*pi;
                dir = [cos(theta);sin(theta)];

                for i = 1:size(dir,2)
                    di = dir(:,i);

                    t = 1/sqrt(di'*obj.P*di);
                    xcur = obj.q + t*di;       
                    x = [x, xcur];     
                end
                center = obj.q;
            
            elseif nargin > 1

                checkBvalidity(obj,B,2);

                if nargin < 3
                    flag = 2;
                end

                switch(flag)
                    case 1
                        e_to_use = obj.sectEll(B);
                    otherwise
                        e_to_use = obj.projEll(B);
                end

                [x, center] = e_to_use.get2DBoundary();

            else
                error('Unable to calculate ellipsoid boundary');
            end
        end


        function EP = projEll(obj,B)
            % This method computes the projection of the ellipsoid onto a 
            % subspace, specified by orthogonal basis vectors B.
            checkBvalidity(obj,B);

            % Computing the projection
            EP = ell(inv(B*inv(obj.P)*B'),B*obj.q);
        end

        function EP = sectEll(obj,B)
            % This method computes the sectioning of the ellipsoid onto a 
            % subspace, specified by orthogonal basis vectors B.
            checkBvalidity(obj,B);

            % Computing the section
            EP = ell(B*obj.P*B',B*obj.q);
        end

        function plotEll(obj,B,flag)
            % This method plots the ellipsoid in 2-D space

            if obj.n == 2

                yalmip('clear');
                x = sdpvar(2,1);
    
                con = [   1     (x-obj.q)'; 
                      (x-obj.q) inv(obj.P)] >= 0;
    
                plot(con);
                hold on;
                plot(obj.q(1),obj.q(2),'o'); 

            elseif nargin > 1
                
                checkBvalidity(obj,B,2);

                if nargin < 3
                    flag = 2;
                end

                switch(flag)
                    case 1
                        e_to_use = obj.sectEll(B);
                    otherwise
                        e_to_use = obj.projEll(B);
                end

                e_to_use.plotEll();

            else

                error('This method plots 2D ellipsoids');
            end
        end

        function plotEllBoundary(obj,B,style,flag)
            % This method plots the boundary of the ellipsoid in 2-D space

            if obj.n == 2

                [x, center] = get2DBoundary(obj);

                plot(center(1),center(2),'*','Color','r');
                hold on;
                
                if nargin <= 2
                    style = 'm';
                end
                
                plot(x(1,:),x(2,:),style);

            elseif nargin > 1

                checkBvalidity(obj,B,2);

                if nargin < 4
                    flag = 2;
                end

                switch(flag)
                    case 1
                        e_to_use = obj.sectEll(B);
                    otherwise
                        e_to_use = obj.projEll(B);
                end

                if nargin > 2
                    e_to_use.plotEllBoundary(NaN,style);
                else
                    e_to_use.plotEllBoundary();
                end
                
            else
                error('This method plots 2D ellipsoids');
            end

        end


    end



end

