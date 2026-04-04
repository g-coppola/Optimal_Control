function [A_cell,B_cell,xeq_used] = build_pol(param,f,xgrid,x,u)

g = param.g;
m = param.M;
n = param.n;
dt = param.dt;

A_tilde = jacobian(f,x);
B_tilde = jacobian(f,u);

A_eq = matlabFunction(A_tilde,'Vars',{x,u});
B_eq = matlabFunction(B_tilde,'Vars',{x,u});

A_list = {};
B_list = {};
xeq_used = [];

tol = 1e-6; 

for i = 1:length(xgrid)
    x1eq = xgrid(i);
    xt = [x1eq, 0, 0, 0]';
    ut = m*g*x1eq;

    At = A_eq(xt,ut);
    Bt = B_eq(xt,ut);

    sys = ss(At,Bt,eye(n),0);
    sysd = c2d(sys, dt, 'zoh');
    Ad = sysd.A;
    Bd = sysd.B;

    if isempty(A_list)
        A_list{end+1} = Ad;
        B_list{end+1} = Bd;
        xeq_used(end+1) = x1eq;
        continue
    end

    nVert = length(A_list);
    A_mat = cell2mat(cellfun(@(x) x(:), A_list, 'UniformOutput', false));
    B_mat = cell2mat(cellfun(@(x) x(:), B_list, 'UniformOutput', false));

    fvec = Ad(:);
    H = [];
    f = [];
    
    Aeq = [A_mat; ones(1,nVert)];
    beq = [fvec; 1];
    lb = zeros(nVert,1);
    ub = [];

    opts = optimoptions('linprog','Display','off');
    [lambda,~,exitflag] = linprog(zeros(nVert,1),[],[],Aeq,beq,lb,ub,opts);

    if exitflag ~= 1 
        A_list{end+1} = Ad;
        B_list{end+1} = Bd;
        xeq_used(end+1) = x1eq;
    end
end

A_cell = A_list;
B_cell = B_list;
end