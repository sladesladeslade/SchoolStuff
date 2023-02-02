clear all; clc()

theta = 0:1:360;

sigman = (22.5)/(.015*.01)*(cos(theta*pi/180)).^2;

taun = [(22.5*cos(theta*pi/180)-22.5*(cos(theta*pi/180)).^3)/(.01*.015) 0 ((-22.5*(cos(theta*pi/180)).^2).*sin(theta*pi/180))/(.01*.015)];

plot(theta, sigman, linewidth=2)
%plot(theta, norm(taun), linewidth=2)