clear all; clc()

theta = 0:1:360;

sigman = (22.5)/(.015*.01)*(cos(theta*pi/180)).^2;

taun = [];

plot(theta, sigman, linewidth=2)
%plot(theta, norm(taun), linewidth=2)