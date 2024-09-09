%% Initialize
clear; clc;

%% Define variables and constants
tspan = [0,9];
y0 = 10;

%% Solve non-stiff differential equations, medium order method.
[t, y] = ode45(@dstate, tspan, y0);

%% Plot results
figure;
clf;
plot(t, y);
grid;
xlabel('t (s)');
ylabel('y (m)');

%% Define derivative function
function dydt = dstate(t,y)
    %% Define constants
    alpha = 2;
    gamma = 0.0001;

    %% Define derivatives
    dydt = alpha*y - gamma*y^2;
end