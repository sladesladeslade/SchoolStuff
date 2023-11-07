clear; clc; close all;

syms n x a b v E C11 C12 C22 C33 P h

% B matrix
B = [-(1-n)/a 0 (1-n)/a 0 (1+n)/a 0 -(1+n)/a 0;
    0 -(1-x)/b 0 -(1+x)/b 0 (1+x)/b 0 (1-x)/b;
    -(1-x)/b -(1-n)/a -(1+x)/b (1-n)/a (1+x)/b (1+n)/a (1-x)/b -(1+n)/a];

% C matrix
C = (E*(1-v))/((1+v)*(1-v))*[1 v/(1-v) 0;
                            v/(1-v) 1 0;
                            0 0 (1-2*v)/(2*(1-v))];

% ke matrix
ke = (a*b*h/16)*int(int(transpose(B)*C*B, x, [-1, 1]), n, [-1, 1]);

% fe matrix
fe = [0; 0; 0; -P; 0; 0; 0; 0];

% calc disp at only node 2
disp2 = inv(ke(3:4, 3:4))*fe(3:4, 1);
fprintf("x disp @ node 2:")
simplify(disp2(1))
fprintf("y disp @ node 2:")
simplify(disp2(2))