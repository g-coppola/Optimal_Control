function [A_cell,B_cell] = build_pol(param,f,xeq_vals,x,u)
g = param.g;
m = param.M;

n = param.n;
dt = param.dt;

A_tilde = jacobian(f,x);
B_tilde = jacobian(f,u);

A_eq = matlabFunction(A_tilde,'Vars',{x,u});
B_eq = matlabFunction(B_tilde,'Vars',{x,u});

nVert = length(xeq_vals);

A_cell = cell(1,nVert);
B_cell = cell(1,nVert);

for i = 1:nVert
    x1eq = xeq_vals{i};
    xt = [x1eq, 0, 0, 0]';
    ut = m*g*x1eq;
    At = A_eq(xt,ut);
    Bt = B_eq(xt,ut);

    sys = ss(At,Bt,eye(n),0);
    sysd = c2d(sys, dt, 'zoh');

    A_cell{i} = sysd.A;
    B_cell{i} = sysd.B; 
end

end

