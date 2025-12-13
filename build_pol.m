function [A_cell,B_cell] = build_pol(param,f,x1_vals,x,u)
g = param.g;
m = param.M;
n = param.n;
dt = param.dt;

A_tilde = jacobian(f,x);
B_tilde = jacobian(f,u);

A_eq = matlabFunction(A_tilde,'Vars',{x,u});
B_eq = matlabFunction(B_tilde,'Vars',{x,u});

nVert = length(x1_vals);

A_cell = cell(1,nVert);
B_cell = cell(1,nVert);

for i = 1:nVert
    xt = [x1_vals{i}, 0, 0, 0]';
    ut = m*g*x1_vals{i}*cos(0);
    At = A_eq(xt,ut);
    Bt = B_eq(xt,ut);

    sys = ss(At,Bt,eye(n),0);
    sysd = c2d(sys, dt, 'zoh');

    A_cell{i} = sysd.A;
    B_cell{i} = sysd.B; 
end

end

