clc; clear all;

theta = 0.7454;
L3 = 1.769;
sigma = 6.24e-4;
sigmap = 6.786e-4;
I = 5.235e-7;
Ip = 4.405e-7;
q = 700;
a = 1.3;
b = 2.7;
L = 4;
E = 73.1e9;
Ep = 68.9e9;

A = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
    a 0 0 0 0 0 0 0 0 -b 0 0 0 0 0 0 0 0;
    0 (a^3)/6 (a^3)/2 0 0 0 0 -(a^3)/6 -(b*(a^2))/2 0 -(b^2)*a -(b^3) 0 0 0 0 0 0;
    (a*cos(theta))/(E*sigma) ((a^3)*sin(theta))/(6*E*I) ((a^3)*sin(theta))/(2*E*I) 0 0 0 0 0 0 0 0 0 -(L3)/(Ep*sigmap) 0 0 0 0 0;
    (a*sin(theta))/(E*sigma) ((a^3)*cos(theta))/(6*E*I) ((a^3)*cos(theta))/(2*E*I) 0 0 0 0 0 0 0 0 0 0 0 0 0 -(L3^3)/(Ep*Ip) 0;
    0 (a^2)/2 a^2 0 0 0 0 -(a^2)/2 -b*a 0 -(b^2) 0 0 0 0 0 0 0;
    -cos(theta) sin(theta) 0 0 0 0 0 -sin(theta) 0 0 0 0 -1 0 0 0 0 0;
    sin(theta) cos(theta) 0 0 0 0 0 -cos(theta) 0 0 0 0 0 0 0 0 0 0;
    0 a a 0 0 0 0 -a -b 0 0 0 0 0 0 0 0 0];
B = [0; 0; 0; 0; 0; 0; 0; 0; -q*L; (q*L^2)/(2*b); 0; 0; -(q*(a^4)*sin(theta))/(24*E*I); -(q*(a^4)*cos(theta))/(24*E*I); 0; 0; 0; 0];

% yikes
Cs = A\B;

La = 0:0.001:a;
Lb = a:0.001:L;
Lt = 0:0.001:L3;
%% Shear
V1 = q*La + Cs(2);
V2 = q*Lb + Cs(8);
plot(La, V1, Lb, V2);

%% Moment
M1 = .5*q*La.^2+Cs(2)*La+Cs(3)*a;
M2 = .5*q*Lb.^2+Cs(8)*Lb+Cs(9)*b;
plot(La, M1, Lb, M2);

%% Rotation
p1 = (1/(E*I))*((1/6)*q*La.^3 + 0.5*Cs(2)*La.^2 + Cs(3)*a*La + Cs(5)*a^2);
p2 = (1/(E*I))*((1/6)*q*Lb.^3 + 0.5*Cs(8)*Lb.^2 + Cs(9)*b*Lb + Cs(11)*b^2);
figure(1);plot(La, p1, Lb, p2);
p3 = (1/(Ep*Ip))*(0.5*Cs(14)*Lt.^2 + Cs(15)*L3*Lt + Cs(17)*L3^2);
%figure(2);plot(Lt, p3);

%% spar vertical displacement
w1 = (1/(E*I))*((1/24)*q*La.^4 + (1/6)*Cs(2)*La.^3 + 0.5*Cs(3)*a*La.^2 + Cs(5)*a^2*La + Cs(6)*a^3);
w2 = (1/(E*I))*((1/24)*q*Lb.^4 + (1/6)*Cs(8)*Lb.^3 + 0.5*Cs(9)*b*Lb.^2 + Cs(11)*b^2*Lb + Cs(12)*b^3);
plot(La, w1, Lb, w2);