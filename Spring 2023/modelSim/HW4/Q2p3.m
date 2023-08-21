clc(); clear all

syms s

A = [-5 3; 1 -4];
B = [0; 5];
C = [1 0; 0 1];
D = [0; 0];

X = C*inv(s*eye(2) - A)*B + D;
X_exp = collect(X, s)