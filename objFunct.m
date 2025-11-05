function [J] = objFunct(Z,h,Q,R,P,param)
    J = 0;
    n = param.n;
    m = param.m;
    
    num_X_vars = n * (h + 1);
    X_vec = Z(1:num_X_vars);
    U_vec = Z(num_X_vars + 1 : end);
    
 
    X = reshape(X_vec, n, h + 1);
    U = reshape(U_vec, m, h);
 
    for k = 1:h
        u = U(:,k); 
        x = X(:,k);
        J = J + x' * Q * x + u' * R * u;
    end

    % Costo terminale
    xf = X(:,h+1); 
    J = J + xf' * P * xf;
end