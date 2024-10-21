%% Initialize
clear; clc;

%% Define variables and constants
mu = 398600;        % km^3/s^2

%% Initial conditions
R = [2500 16000 4000];     % km
V = [-3 -1 5];             % km/s

%% run matlab script from appendix
coe = coe_from_sv(R, V, mu);

% outputs
% coe = [h e RA incl w TA a]
fprintf("e=%g", coe(2));
fprintf("\nh=%g km^2/s", coe(1));
fprintf("\ni=%g deg", rad2deg(coe(4)));
fprintf("\nOmega=%g deg", rad2deg(coe(3)));
fprintf("\nw=%g deg", rad2deg(coe(5)));
fprintf("\ntheta=%g deg", rad2deg(coe(6)));