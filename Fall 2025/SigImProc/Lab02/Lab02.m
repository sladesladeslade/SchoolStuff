%%
% Slade Brooks
% brooksl@mail.uc.edu
% 09.03.25
% BME6013C
% Lab 02

clear variables
close all

%% Part 1
% define Gaussian function coefs
% g = A*exp(-(x^2 + y^2)/(2*sigma^2)
A = 255;
sigma = 10;

% define x and y range (abs <= 5sigma)
step = 0.1;
x = -5*sigma:step:5*sigma;
y = x;

% meshgrid data points
[X, Y] = meshgrid(x, y);

% get vals for gaussian function
g = A*exp(-(X.^2 + Y.^2)/(2*sigma^2));

% plotting
figure()
imagesc(x, y, g)
axis equal;
c = colorbar(); c.Label.String = "Amplitude (cm)";
title("Gaussian Function"); xlabel("x (cm)"); ylabel("y (cm)");
xticks(x(1):10:x(end)); yticks(y(1):10:y(end));
xlim([x(1) x(end)]); ylim([y(1) y(end)]);

%% Part 2
% discretized fxn
gdisc = round(g);

% calculate noise
n = g - gdisc;

% plotting
figure()
imagesc(x, y, n)
axis equal;
c = colorbar(); c.Label.String = "Amplitude (cm)";
title("Noise of Discretized Gaussian Function"); xlabel("x (cm)"); ylabel("y (cm)");
xticks(x(1):10:x(end)); yticks(y(1):10:y(end));
xlim([x(1) x(end)]); ylim([y(1) y(end)]);

%% Part 3
% calc and disp mean and rms of both
g_mean = mean(g, "All"); disp("g_mean: " + g_mean);
n_mean = mean(n, "All"); disp("n_mean: " + n_mean);
g_rms = rms(g, "All"); disp("g_rms: " + g_rms);
n_rms = rms(n, "All"); disp("n_rms: " + n_rms);

%% Part 4
% calculate SNR
SNR_dB = 20*log10(g_rms/n_rms);
disp("SNR_dB: " + SNR_dB)