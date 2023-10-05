clc; clear all; close all

% Problem 1
% Create the stiffness matrix K
K = 10^8 * [
    2 0 -2 0 0 0;
    0 0 0 0 0 0;
    -2 0 (2+1/sqrt(2)) (-1/sqrt(2)) (-1/sqrt(2)) (1/sqrt(2));
    0 0 (-1/sqrt(2)) (1/sqrt(2)) (1/sqrt(2)) (-1/sqrt(2));
    0 0 (-1/sqrt(2)) (1/sqrt(2)) (1/sqrt(2)) (-1/sqrt(2));
    0 0 (1/sqrt(2)) (-1/sqrt(2)) (-1/sqrt(2)) (1/sqrt(2))
];
% Create the force vector F
F = [
    0;
    -100;
];
% Solve for D for the 3rd and 4th rows using the equation K * D = F
D = K(3:4, 3:4) \ F;

% solve for F1
F1 = (0.01 * 2E10)/1 * [-1 0 1 0] * [0; 0; -5E-7; -1.9142E-6];
F2 = (0.01 * 2E10)/sqrt(2) * [-1/sqrt(2), 1/sqrt(2), 1/sqrt(2), -1/sqrt(2)] * [-5E-7; -1.9142E-6; 0; 0];

% create new D
De = [0; 0; D(1); D(2); 0; 0];
% solve for each F from F
Fe = K*De;