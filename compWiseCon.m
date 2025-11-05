function [F,Qe] = compWiseCon(A,B,C,x0,y_max,u_max)
[n,m]=size(B{1});
ny = size(C,1);
Q1=eye(n);
R=eye(m);
yalmip('clear')
Q=sdpvar(n);
gamma=sdpvar;
Y=sdpvar(m,n,'full');
con=[Q>=0, [1 x0';x0 Q]>=0];
for j = 1:length(A)
    Aj = A{j};
    Bj = B{j};
    con=[con, [Q Q*Aj'+Y'*Bj' Q*sqrt(Q1) Y'*sqrt(R);
                Aj*Q+Bj*Y Q zeros(n) zeros(n,m);
                sqrt(Q1)*Q zeros(n) gamma*eye(n) zeros(n,m);
                sqrt(R)*Y zeros(m,n) zeros(m,n) gamma*eye(m)]>=0];
    for l = 1:ny
        Cl = C(l,:);
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

