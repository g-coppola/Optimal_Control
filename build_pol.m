function [A_cell,B_cell] = build_pol(param,dt,L,th)
s_max = 1;
s_min = sin(th)/th;

c_max = cos(0);
c_min = cos(th);

I = param.I;
g = param.g;
m = param.M;
M = param.MM;
n = param.n;

d_max = 1/I;
d_min = 1/(I+M*L^2);

num_params = 3;
num_vert = 2^num_params;

A_cell = cell(1, num_vert);
B_cell = cell(1, num_vert);



for i = 0:(num_vert - 1)
    bin_comb = dec2bin(i, num_params);
    
    if bin_comb(3) == '0', phi1_val = s_min; else, phi1_val = s_max; end
    if bin_comb(2) == '0', phi2_val = c_min; else, phi2_val = c_max; end
    if bin_comb(1) == '0', phi3_val = d_min; else, phi3_val = d_max; end
    
    A_temp = zeros(4, 4);
    A_temp(1, 2) = 1;
    A_temp(2, 3) = -(m*g/M) * phi1_val;
    A_temp(3, 4) = 1;
    A_temp(4, 1) = -m*g * phi2_val * phi3_val;
    
    B_temp = zeros(4, 1);
    B_temp(4, 1) = phi3_val;
    
    sysc = ss(A_temp,B_temp,eye(n),0);
    ssd = c2d(sysc, dt, 'zoh');

    A_cell{i+1} = ssd.A;
    B_cell{i+1} = ssd.B;

    
end
end

