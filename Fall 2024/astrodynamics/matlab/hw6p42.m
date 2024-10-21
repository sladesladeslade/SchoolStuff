%% Initialize
clear; clc;

%% Define variables and constants
Re = 6378;          % km
global mu
mu = 398600;        % km^3/s^2

%% given initial conditions
h = 500;            % km
a0 = 300;           % deg
d0 = -20;           % deg

%% determine initial location and velocity
R0 = (6378 + h)*[cosd(d0)*cosd(a0) cosd(d0)*sind(a0) sind(d0)];     % km
V0 = [0 0 10];                                                      % km/s
t = 1800;                                                           % s

%% use algorithm from appendix again
[R, V] = rv_from_r0v0(R0, V0, t);

%% use algortihm 4.1
% get these vars
r = norm(R);        % km
l = R(1)/r;
m = R(2)/r;
n = R(3)/r;

% calc and return del
d = asind(n);
fprintf("del: %g deg", d);

% check value of m
fprintf("\nm: %g", m);

% b/c m>0
a = acosd(l/cosd(d));
fprintf("\nalpha: %g deg", a);