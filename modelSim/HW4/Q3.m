clc; clear all

% make matrices for vals
A = [0 1; -4/5 -7/5];
B = [0; 1/5];
C = [1 0];
D = [0];

% create sys
sys = ss(A, B, C, D);

% set initial conditions
x0 = [1 0];

% plot response
figure(1)
initial(sys, x0);

clc; clear all

% call ode45
[t, y] = ode45(@xdot_prob, [0, 10], [1, 0]);
figure(2)
plot(t, y)
xlabel("Time (s)"); ylabel("y"); legend("$\dot{x_{1}}$", "$\dot{x_{2}}$", "interpreter", "latex")

% function for ODE
function xdot = xdot_prob(t, x)
    xdot = [x(2); -4/5*x(1)-7/5*x(2)];
end