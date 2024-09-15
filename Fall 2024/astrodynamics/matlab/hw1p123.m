%% Initialize
clear; clc;

%% Define variables and constants
tspan = [0,100];
x0 = 0;
y0 = 1;
z0 = 0;
ics = [x0; y0; z0];

%% Solve non-stiff differential equations, medium order method.
[t, y] = ode45(@dstate, tspan, ics);

%% Plot results
figure;
clf;
plot3(y(:,1), y(:,2), y(:,3));
grid;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

%% Define derivative function
function ddt = dstate(t, yi)
    %% Define constants
    sigma = 10;
    beta = 8/3;
    rho = 28;
    
    %% get state from inputs
    x = yi(1);
    y = yi(2);
    z = yi(3);

    %% Define derivatives
    dxdt = sigma*(y - x);
    dydt = x*(rho - z) - y;
    dzdt = x*y - beta*z;
    
    %% create output vector
    ddt = [dxdt; dydt; dzdt];
end