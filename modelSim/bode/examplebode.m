clear; clc

%% Parameters

R = 1; C = 0.01; L = 0.01;

%% Low Pass Filter

L_RC = tf(1, [R*C 1]);
L_LR = tf(R, [L R]);
figure(1); bode(L_RC); grid on
figure(2); bode(L_LR); grid on

%% High Pass Filter