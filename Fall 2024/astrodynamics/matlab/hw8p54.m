%% Initialize
clear; clc;

%% Problem 5.4
% Define variables and constants
global mu
mu = 398600;        % km^3/s^2

% Initial conditions
r1 = [3600 4600 3600];          % km
r2 = [-5500 6240 -5200];        % km
dt = 30*60;                     % s
string = "pro";

% use lambert function (alogrithm 5.2) from book to get v at each position
[v1, v2] = lambert(r1, r2, dt, string);

%% calculate epsilon
% get mangitude of v and r
v = norm(v1);
r = norm(r1);

% calculate and return epsilon
e = v^2/2 - mu/r;
fprintf("epsilon=%g km^2/s^2", e)

%% Problem 5.5
% calculate h
h1 = cross(r1, v1);
h = norm(h1);

% calc e
e = sqrt(1 + h^2/mu^2*(v^2 - 2*mu/r));

% then get r, z, and i
r = h^2/mu*1/(1 + e);
zp = r - 6378;
i = acosd(h1(3)/h);

% print results
fprintf("\nz perigee=%g km", zp);
fprintf("\ni=%g deg", i);