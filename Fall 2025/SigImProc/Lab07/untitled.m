% BME5013C/6113C: Biomedical Signal and Image Processing
% Fall 2025
% Fourier Transform Window Example Set #02
%
% Author : T. Douglas Mast, PhD
% Andres Llico, PhD (edits)
% Date : Oct 06, 2025

clear variables;
close all;

% Define a spatial grid for all examples
xmax = 4; dx = 0.01; % spatial range and step (m)
ymax = xmax; dy = dx; % vertical dimensions matching horizontal
xvec = -xmax:dx:xmax;
yvec = -ymax:dy:ymax;
[x,y] = meshgrid(xvec,yvec);

% Define a spatial-frequency grid as well
umax = 2; du = 0.004;
vmax = umax; dv = du;
uvec = -umax:du:umax;
vvec = -vmax:dv:vmax;
[u,v] = meshgrid(uvec,vvec);

%% Example 1: plane wave, spatial frequency nu0, direction phi
nu0 = 1.5; phi = +pi/3;
f = cos(2*pi*nu0*cos(phi)*x + 2*pi*nu0*sin(phi)*y);

% Approximately represent delta functions in F.T. by Gaussian spikes
sigma = 0.025;
F = ( exp(-(u-nu0*cos(phi)).^2/(2*sigma.^2) ...
- (v-nu0*sin(phi)).^2/(2*sigma.^2)) ...
+ exp(-(u+nu0*cos(phi)).^2/(2*sigma.^2) ...
- (v+nu0*sin(phi)).^2/(2*sigma.^2)) ) / (4*pi*sigma^2) ;
figure(1)
set(gcf, ...
"Units", "inches", ...
"Position", [4 1 8 3]);
subplot(1,2,1);
imagesc(xvec, yvec, f, [-1 1]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("x (m)", ...
"FontSize", 14);
ylabel("y (m)", ...
"FontSize", 14);
title("Plane wave f(x,y)", ...
"FontSize", 16);
subplot(1,2,2);
imagesc(uvec, vvec, F, [0 max(max(F))]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("u (m^{-1})", ...
"FontSize", 14);
ylabel("v (m^{-1})", ...
"FontSize", 14);
title("Delta functions F(u,v)", ...
"FontSize", 16);

%% Example 2: Box function (F.T. sinc function)
Lx = 1.; Ly = 4.5;
f = ((abs(x)<=Lx/2)&(abs(y)<=Ly/2)) / (Lx*Ly);
F = sinc(u*Lx).*sinc(v*Ly);
figure(2);
set(gcf, ...
"Units", "inches", ...
"Position", [4 1 8 3]);
subplot(1,2,1);
imagesc(xvec,yvec,f,[0 max(max(f))]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
set(gca,'YDir','Normal');
xlabel("x (m)", ...
"FontSize", 14);
ylabel("y (m)", ...
"FontSize", 14);
title("Box function f(x,y)", ...
"FontSize", 16);
subplot(1,2,2);
imagesc(uvec,vvec,F,[-1 1]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
set(gca,'YDir','Normal');
xlabel("u (m^{-1})", ...
"FontSize", 14);
ylabel("v (m^{-1})", ...
"FontSize", 14);
title("Sinc function F(u,v)", ...
"FontSize", 16);

%% Example 3: Gaussian (F.T. also Gaussian)
sigmax = 1; sigmay = .3;
f = exp(-x.^2/(2*sigmax^2)).*exp(-y.^2/(2*sigmay^2))/(2*pi*sigmax*sigmay);
F = exp(-2*pi^2*(u.^2*sigmax^2 + v.^2*sigmay^2));
figure(3);
set(gcf, ...
"Units", "inches", ...
"Position", [4 1 8 3]);
subplot(1,2,1);
imagesc(xvec,yvec,f,[0 max(max(f))]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("x (m)", ...
"FontSize", 14);
ylabel("y (m)", ...
"FontSize", 14);
title("Gaussian f(x,y)", ...
"FontSize", 16);
subplot(1,2,2);
imagesc(uvec, vvec, F, [0 1]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("u (m^{-1})", ...
"FontSize", 14);
ylabel("v (m^{-1})", ...
"FontSize", 14);
title("Gaussian F(u,v)", ...
"FontSize", 16);

%% Example 4: Circ (F.T. Jinc function)
r0 = 1; r = sqrt(x.^2+y.^2);
f = (r <= r0)/(pi*r0^2);
nu = sqrt(u.^2+v.^2);
F = besselj(1, 2*pi*nu*r0)./(pi*nu*r0);
F(nu == 0) = 1;
figure(4);
set(gcf, ...
"Units", "inches", ...
"Position", [4 1 8 3]);
subplot(1,2,1);
imagesc(xvec, yvec, f, [0 max(max(f))]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("x (m)", ...
"FontSize", 14);
ylabel("y (m)", ...
"FontSize", 14);
title("Circ f(x,y)", ...
"FontSize", 16);
subplot(1,2,2);
imagesc(uvec, vvec, F, [-1 1]); colorbar;
colormap(gray); axis image;
set(gca,'YDir','Normal');
xlabel("u (m^{-1})", ...
"FontSize", 14);
ylabel("v (m^{-1})", ...
"FontSize", 14);
title("Jinc F(u,v)", ...
"FontSize", 16);