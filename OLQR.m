function [PP,FF]=OLQR(A,B,Q,R,S,N)

PP = cell(1,N);
PP{N} = S;
FF = cell(1,N);

    for k = N-1:-1:1
        FF{k+1} = inv(R + B'*PP{k+1}*B)*B'*PP{k+1}*A;
        PP{k} = A'*PP{k+1}*A + Q - A'*PP{k+1}*B*inv(R+B'*PP{k+1}*B)*B'*PP{k+1}*A;
    end

FF{1} = inv(R + B'*PP{2}*B)*B'*PP{2}*A;
end