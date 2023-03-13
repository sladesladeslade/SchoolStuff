clc; clear all
%% 
% call ode45
[t, y] = ode45(@HW63, [0, 300], 0);
plot(t, y)
title("ODE45 - HW6.3");
xlabel("Time (s)"); ylabel("V (m/s)");

% function for ODE
function vdot = HW63(t, v)
    I = 1.1;
    m = 50;
    R = 0.3;
    T = 500;
    ct = 0.2;
    g = 9.81;
    vdot = (T*R - m*g*R^2 - ct*v)/(I + m*R^2);
end

%%