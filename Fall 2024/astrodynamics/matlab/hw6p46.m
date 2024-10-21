%% Initialize
clear; clc;

%% given vals
R = [-6000 -1000 -5000];        % km
e = [0.4 0.5 0.6];

%% determine magnitudes
r = norm(R);
em = norm(e);

%% calculate true anomaly
% approaching perigee so:
theta = 360 - acosd(dot(R, e)/(em*r));
fprintf("true anomaly: %g deg", theta);