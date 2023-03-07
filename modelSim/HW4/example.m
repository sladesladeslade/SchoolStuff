[t, x] = ode45(@fxn1, [0, 7], 0);
figure(1); plot(t, x)

[t, y] = ode45(@fxn2, [0, 3.5], [0;0]);
figure(2); plot(t, y); legend("dumb", "dumber")

function xdot = fxn1(t, x)
    f = sin(t);
    xdot = -7/5*x + 3*f;
end

function ydot = fxn2(t, y)
    f = 1;
    ydot = [ y(2); -21*y(1)-4*y(2)+4*f];
end