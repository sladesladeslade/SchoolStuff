%% Initialize
clear; clc;

%% Define variables and constants
tspan = [0 24*3600];
mu = 398600;
R = 6378;
r0 = [6600 0 0];
v0 = [0 12 0];
ics = [r0'; v0'];

%% Solve non-stiff differential equations, medium order method.
[t, y] = ode45(@dstate, tspan, ics);

%% print results
fprintf("dist. @ 24hr: %.2f km \n", norm([y(end,1) y(end,2) y(end,3)]));
fprintf("speed @ 24hr: %.2f km/s", norm([y(end,4) y(end,5) y(end,6)]));

%% Plot results
%figure;
%clf;
%plot3(y(:,1), y(:,2), y(:,3));
%grid;
%xlabel('x (m)');
%ylabel('y (m)');
%zlabel('z (m)');

%% Define derivative function
function ddt = dstate(t, yi)  
    %% define constants
    mu = 398600;

    %% get state from inputs
    x = yi(1);
    y = yi(2);
    z = yi(3);
    vx = yi(4);
    vy = yi(5);
    vz = yi(6);
    r = norm([x y z]);

    %% Define derivatives
    ax = -mu*x/r^3;
    ay = -mu*y/r^3;
    az = -mu*z/r^3;
    
    %% create output vector
    ddt = [vx; vy; vz; ax; ay; az];
end