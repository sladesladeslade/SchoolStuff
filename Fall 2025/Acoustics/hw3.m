%%
% Slade Brooks
% brooksl@mail.uc.edu
% 10.19.25
% MECH6066
% HW 03

clear variables
close all

%% set up vals
% range of frequencies
fs = 0.1e6:0.01e6:10e6;
f = 2e6;

% acoustic properties
Z1 = 1.48e6; c1 = 1481;
Z2 = 17e6; c2 = 6300;
ZL = sqrt(Z1*Z2); cL = sqrt(c1*c2); L = cL/(4*f);
Zs = 3.58e6; cs = 1730; Ls = cs/(4*f);

%% calculate transmission
% w/o matching layer
TIi = 4*Z1*Z2/(Z1 + Z2)^2;

% w/ ideal matching layer
w = (2*pi).*fs;
K = w./cL;
TIii = 4./(2 + (Z2/Z1 + Z1/Z2).*(cos(K*L)).^2 + ((ZL^2)/(Z1*Z2) + (Z1*Z2)/(ZL^2)).*(sin(K*L)).^2);

% w/ practical matching layer
w = (2*pi).*fs;
K = w./cs;
TIiii = 4./(2 + (Z2/Z1 + Z1/Z2).*(cos(K*Ls)).^2 + ((Zs^2)/(Z1*Z2) + (Z1*Z2)/(Zs^2)).*(sin(K*Ls)).^2);

%% plotting
figure
hold on
yline(TIi); plot(fs, TIii); plot(fs, TIiii);
hold off
xlabel("Frequency (Hz)"); xlim([fs(1) fs(end)]);
ylabel("Intensity Transmission Coefficient"); ylim([0 1]);
legend(["w/o matching layer" "w/ ideal matching layer" "w/ practical matching layer"], location="southeast");
grid on