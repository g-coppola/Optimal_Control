function [c,ceq] = nonLinConstraints(Z,xk,A,B,h,P,dt,ft,param)
    c = [];
    ceq = [];

    n = param.n;
    m = param.m;
    l = size(A,2);
    
    num_X_vars = n * (h + 1);
    X_vec = Z(1:num_X_vars);
    U_vec = Z(num_X_vars + 1 : end);
    
    X = reshape(X_vec, n, h + 1);
    U = reshape(U_vec, m, h);
    
    for i = 1:l
        cf = (A{i}*X(:,h)+B{i}+U(:,h))'*P*(A{i}*X(:,h)+B{i}+U(:,h))-1;
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