function [J] = objFunct(Z,h,A,B,Q,R,P,param)
    J = 0;
    n = param.n;
    m = param.m;
    l = size(A,2);
    
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
            c(i) = (A{i}*x+B{i}*u)'*Q*(A{i}*x+B{i}*u)+ u'*R *u;
        end
        J = J + max(c);
    end

    % Terminal Cost
    cf = zeros(1,l);
    for i = 1:l
        cf(i) = (A{i}*X(:,h)+B{i}*U(:,h))'*P*(A{i}*X(:,h)+B{i}*U(:,h));
    end
    J = J + max(cf);
end