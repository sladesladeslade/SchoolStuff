%%
% Slade Brooks
% brooksl@mail.uc.edu
% 09.10.25
% BME6013C
% Lab 03

%clear variables
%close all

%% Part 1
%% Part 1.a
% define time step and length
T = 1;
dt = 0.0001;
t = 0:dt:T;

% define m and n
m = 1;
n = 5;

% calculate the signal
g0 = exp(-j*2*pi*m*t/T).*exp(j*2*pi*n*t/T);

% determine the real, imaginary, and magnitude
g0r = real(g0);
g0im = imag(g0);
g0mag = abs(g0);

% plot data
figure()
hold on
plot(t, g0r, linewidth=1.5); plot(t, g0im, linewidth=1.5); plot(t, g0mag, linewidth=1.5);
hold off
legend({"re\{g0\}","im\{g0\}","|g0|"}); title("Plot of signal components and magnitude")
ylim([-1.5 1.5]); xlim([0 T]); ylabel("Amplitude (cm)"); xlabel("t (s)"); grid on;

% using datatips, measure the signal at a couple peaks (time btn 3 peaks)
ts = 0.7495 - 0.2498;     % measured time diff
npks = 2;                 % n cycles
f_est = npks/ts;
disp("f measured = " + f_est + " Hz") % there is a frequency of 4 Hz visible in the data

%% Part 1.b
% calculate rfg (the numerical approx should be the mean of a discrete signal)
% so rfg = 1/N sum(g(0: N-1))
r_fg_1b = 1/size(t, 2)*sum(g0(1:end - 1));
disp("rfg = " + r_fg_1b);      % r_fg is e-17 so basically 0

%% Part 1.c
% set equal m and n
mc = 5;
nc = 5;

% calculate new values for g0
g0c = exp(-j*2*pi*mc*t/T).*exp(j*2*pi*nc*t/T);

% same as before to determine rfg
r_fg_1c = 1/size(t, 2)*sum(g0c(1:end-1));
disp("rfg c = " + r_fg_1c);   % r_fg is 1 this time as expected

%% Part 2
% load problem file
p_6_5 = load("p_6_5.mat");

% get variable list
whos -file p_6_5

%% Part 2.a
% pull out data array
x = p_6_5.x;
lx = size(x, 2);

% set up new range of Ts
fs = 100;
Ts = 1/fs;
t = linspace(0, Ts*lx, lx);

% plot data
figure()
plot(t, x);
xlim([0 Ts*lx]); ylim([-100 100]);
xlabel("t (s)"); ylabel("Amplitude (cm)"); grid on;
title("Plot of EEG data")

%% Part 2.b
% use xcorr to get r vals and lag vals
[r, lag] = xcorr(x, "coeff");

% plot autocorr
figure()
plot(lag*Ts, r);
xlim([-10 10]); ylim([-0.2 1]);
xlabel("Lag (s)"); ylabel("r"); grid on;
title("Plot of EEG autocorrelation function")

%% Part 2.c
% measure freq from autocorr using datatips
ts = (159 - 21)*Ts;      % measured times lag in samples and convert to seconds
npks = 20;          % number of peaks btn peaks
f = npks/ts;        % calculated frequency
disp("f = " + f + " Hz")

% to verify, check a shorter time span also (same method)
ts2 = (49 - 21)*Ts;
npks2 = 4;
f2 = npks2/ts2;
disp("less accurate f = " + f2 + " Hz")