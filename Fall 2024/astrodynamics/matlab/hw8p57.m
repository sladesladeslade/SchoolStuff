%% Initialize
clear; clc;

%% Define variables and constants
global mu
mu = 398600;        % km^3/s^2

%% Initial conditions
r1 = [5644 -2830 -4170];        % km
r2 = [-2240 7320 -4980];        % km
dt = 20*60;                     % s
string = "pro";

%% use lambert function from book to get v at each position
[v1, v2] = lambert(r1, r2, dt, string);

%% get orbital params from coe from sv script from book
coe = coe_from_sv(r1, v1, mu);

% outputs
% coe = [h e RA incl w TA a]
fprintf("e=%g", coe(2));
fprintf("\nh=%g km^2/s", coe(1));
fprintf("\ni=%g deg", rad2deg(coe(4)));
fprintf("\nOmega=%g deg", rad2deg(coe(3)));
fprintf("\nw=%g deg", rad2deg(coe(5)));
fprintf("\ntheta=%g deg", rad2deg(coe(6)));

%% calculate perigee
rp = coe(1)^2/mu*1/(1+coe(2));       % radius of perigee (km)
zp = rp - 6378;                      % alt of perigee (km)
fprintf("\nz perigee=%g km", zp);