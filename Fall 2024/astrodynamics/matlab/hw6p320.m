%% Initialize
clear; clc;

%% Define variables and constants
global mu
mu = 398600;        % km^3/s^2

%% Initial conditions
R0 = [20000 -105000 -19000];        % km
V0 = [0.9 -3.4 -1.5];               % km/s
t = 2*3600;                         % s

%% do algortihm (from appendix D.16) and output results
[R, V] = rv_from_r0v0(R0, V0, t);
fprintf("Final position: R=(%g, %g, %g) km", R(1), R(2), R(3));
fprintf("\nFinal velocity: V=(%g, %g, %g) km/s", V(1), V(2), V(3));