%%
% Slade Brooks
% brooksl@mail.uc.edu
% 09.24.25
% BME6013C
% Lab 05

clear variables
close all

%% Part 1
% set period T, time array for 2 periods, and amplitude
T = 10;
t = 0:0.01:2*T;
A = 10;

% set first limit of n and corresponding range of vals
N = 3;
ns = 1:2:N;

% create empty array for u vals
us = zeros(N, length(t));

% loop through ns to get each component
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end

% sum all columns to get final set of vals
u = sum(us);

% create a figure and add a subplot for the first N
figure;
subplot(3, 1, 1);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-15 15]); xlim([0 t(end)]);
set(gca, 'ytick', -15:5:15);

% repeat previous steps for next value of N
N = 31;
ns = 1:2:N;
us = zeros(N, length(t));
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end
u = sum(us);

% plotting
subplot(3, 1, 2);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-15 15]); xlim([0 t(end)]);
set(gca, 'ytick', -15:5:15);

% repeat previous steps for next value of N
N = 301;
ns = 1:2:N;
us = zeros(N, length(t));
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end
u = sum(us);

% plotting
subplot(3, 1, 3);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-15 15]); xlim([0 t(end)]);
set(gca, 'ytick', -15:5:15);

% It is clear from N=3 that the square wave is poorly reconstructed. The
% slope is fairly low at the 0 crossing points, not a sharp vertical line.
% The value also oscillates around the correct amplitude, over- and
% undershooting by ~2 cm. N=31 is a better approximation. The slope at the
% 0 crossings is much higher, however, lots of oscillations are visible.
% The value does not overshoot as much as the N=3 case, but the
% oscillations are still visible. Finally, the N=301 case approximates very
% well. The slope at zero crossings appears to be nearly vertical, and the
% value is almost exactly constant at the correct amplitude. However, there
% are very sharp corners where it overshoots and oscillates. These were
% slightly visible in the N=31 case, but are much clearer for N=301.

%% Part 2
% we will just copy the code from the previous section for this, but will
% update ns to go from 1:N in steps of 1 instead of steps of 2

% set first limit of n and corresponding range of vals
N = 3;
ns = 1:1:N;

% create empty array for u vals
us = zeros(N, length(t));

% loop through ns to get each component
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end

% sum all columns to get final set of vals
u = sum(us);

% create a figure and add a subplot for the first N
figure;
subplot(3, 1, 1);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-30 30]); xlim([0 t(end)]);
set(gca, 'ytick', -30:10:30);

% repeat previous steps for next value of N
N = 31;
ns = 1:1:N;
us = zeros(N, length(t));
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end
u = sum(us);

% plotting
subplot(3, 1, 2);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-30 30]); xlim([0 t(end)]);
set(gca, 'ytick', -30:10:30);

% repeat previous steps for next value of N
N = 301;
ns = 1:1:N;
us = zeros(N, length(t));
for i = 1:length(ns)
    us(i,:) = 4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end
u = sum(us);

% plotting
subplot(3, 1, 3);
plot(t, u); grid on;
title("N = " + N); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-30 30]); xlim([0 t(end)]);
set(gca, 'ytick', -30:10:30);

% Including the even integers creates a sawtooth wave instead of a square
% wave. The sawtooth has the same period as our original square wave, but
% its maximum value is double the amplitude of our square wave. We also see
% similar effects based on the value of N, where there are many
% oscillations at low counts of N, which get smoothed out as N increases.
% We also see the same effect of slope at 0 crossing increasing with N.
% Finally, the sharp overshoot at high values of N is also visible for the
% sawtooth. It can also be noted that the sawtooth undershoots the peak
% amplitude at low values of N, but the square wave overshoots the peak at
% low values of N.

%% Part 3
% first we will grab the N definition from before
N = 301;
ns = 1:2:N;

% choose value of n1
n1 = 67;

% set up empty array for storing/summing
us = zeros(N, length(t));
for i = 1:length(ns)
    % add the filtering term to this calculation
    us(i,:) = exp(-ns(i)/n1)*4*A/(ns(i)*pi)*sin(2*pi*ns(i)*t/T);
end
u = sum(us);

% plotting
figure;
plot(t, u); grid on;
title("N = " + N + ", n1 = " + n1); ylabel("Amplitude (cm)"); xlabel("Time (s)");
ylim([-15 15]); xlim([0 t(end)]);
set(gca, 'ytick', -15:5:15);

% The filtered square wave has the same period as the original square wave.
% The amplitude is almost identical, but the filtered square wave does not
% quite hit the full amplitude of the original square wave, it is low by
% ~1%. The biggest change with the filter is the reduction of the sharp
% overshoots at the start and end of each portion where the signal
% amplitude is not crossing 0. The filtering smoothed out the corners of
% the square wave, leaving a gentle slope around to the correct value.
% There are no longer any oscillations in amplitude visible. The slope of
% the 0 crossing appears to be slightly less vertical as the corner
% approaches.