%%
% Slade Brooks
% brooksl@mail.uc.edu
% 08.27.25
% BME6013C
% Lab 01

clear variables
close all

%% file i/o
% load problem file
p_1_2 = load('p_1_2.mat');

% get variable list
whos -file p_1_2
vars = whos('-file', 'p_1_2.mat');

%% sampling signal
% define sampling rate (Hz) and calculate sample period (s)
fs = 10;
Ts = 1/fs;

% determine length of x1 from stored vars
lx1 = vars(1).size(2);

% set up array of time values from 0 to whatever final time based on length
% of x1
t = linspace(0, Ts*lx1, lx1);

% pull out x1 vals from file
x1 = p_1_2.x1;

%% plotting
figure()
plot(t, x1);
xlabel('Time (s)'); ylabel('Signal Amplitude (units)');
grid('on');

%% analysis
% manually measured sample period using data tips
T = 1.845 - 0.7175;

% calculate frequency from period
f = 1/T;
f

%% signals 2 and 3
% get lengths of the vars
lx2 = vars(2).size(2);
lx3 = vars(3).size(2);

% set up time arrays for both
t2 = linspace(0, Ts*lx2, lx2);
t3 = linspace(0, Ts*lx3, lx3);

% load data
x2 = p_1_2.x2;
x3 = p_1_2.x3;

%% plotting
figure()
plot(t2, x2);
xlabel('Time (s)'); ylabel('Signal Amplitude (units)');
grid('on');

figure()
plot(t3, x3);
xlabel('Time (s)'); ylabel('Signal Amplitude (units)');
grid('on');

%% analysis
% manually measure points to determine slope
slope2 = (32.6453 - 7.73446)/(6.48 - 1.4175);
slope3 = (720.751 - 160.946)/(6.885 - 1.31625);
slope2
slope3


















