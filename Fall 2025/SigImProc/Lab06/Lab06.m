%%
% Slade Brooks
% brooksl@mail.uc.edu
% 10.01.25
% BME6013C
% Lab 06

clear variables
close all

%% Part 1
% define center frequency (Hz)
f0 = 100;

% determine period of one cycle
T0 = 1/f0;

% get time for 50 cycles
tf = T0*50;

% set number of samples and calculate time step
N = 10000;
dt = tf/N;

% create time vector centered around 0 with given time step
t = -tf/2:dt:tf/2;

%% Part 2
% create figure
figure;

% define T value
T1 = 5/f0;

% determine heaviside fxn values
hin = T1/2 - abs(t);
H = zeros(size(hin));
H(hin >= 0) = 1;

% calculate u(t) for defined T
u1 = H.*exp(j*2*pi*f0.*t);

% get real component and magnitude of u
realu1 = real(u1);
magu1 = abs(u1);

% plot
subplot(3, 1, 1);
hold on;
plot(t, realu1); plot(t, magu1); grid on;
hold off;
title("T = " + T1 + "s (5/f0)"); ylabel("Amplitude (cm)"); xlabel("Time (s)");
xlim([t(1) t(end)]); ylim([-1.5 1.5]);
set(gca, 'ytick', -1.5:0.5:1.5);
legend("Re\{u(t)\}", "|u(t)|");

% now copy that section and plot again with a new choice for T
% define T value
T2 = 15/f0;

% determine heaviside fxn values
hin = T2/2 - abs(t);
H = zeros(size(hin));
H(hin >= 0) = 1;

% calculate u(t) for defined T
u2 = H.*exp(j*2*pi*f0.*t);

% get real component and magnitude of u
realu2 = real(u2);
magu2 = abs(u2);

% plot
subplot(3, 1, 2);
hold on;
plot(t, realu2); plot(t, magu2); grid on;
hold off;
title("T = " + T2 + "s (15/f0)"); ylabel("Amplitude (cm)"); xlabel("Time (s)");
xlim([t(1) t(end)]); ylim([-1.5 1.5]);
set(gca, 'ytick', -1.5:0.5:1.5);
legend("Re\{u(t)\}", "|u(t)|");

% and one final time
% define T value
T3 = 35/f0;

% determine heaviside fxn values
hin = T3/2 - abs(t);
H = zeros(size(hin));
H(hin >= 0) = 1;

% calculate u(t) for defined T
u3 = H.*exp(j*2*pi*f0.*t);

% get real component and magnitude of u
realu3 = real(u3);
magu3 = abs(u3);

% plot
subplot(3, 1, 3);
hold on;
plot(t, realu3); plot(t, magu3); grid on;
hold off;
title("T = " + T3 + "s (35/f0)"); ylabel("Amplitude (cm)"); xlabel("Time (s)");
xlim([t(1) t(end)]); ylim([-1.5 1.5]);
set(gca, 'ytick', -1.5:0.5:1.5);
legend("Re\{u(t)\}", "|u(t)|");

%% Part 3
% define f array, centered around f0 with some width on either side
f = (f0 - 200):0.1:(f0 + 200);

% calculate the fourier transform for each value of T
Uf1 = T1*sinc((f - f0)*T1);
Uf2 = T2*sinc((f - f0)*T2);
Uf3 = T3*sinc((f - f0)*T3);

% calculate in dB for all with a reference value (using 20log10 instead of
% squaring both terms in the log10)
ref = 0.0001;
db1 = 20*log10(abs(Uf1)/ref);
db2 = 20*log10(abs(Uf2)/ref);
db3 = 20*log10(abs(Uf3)/ref);

% create figure
figure;

% plot each subplot
subplot(3, 1, 1);
plot(f, db1); grid on;
title("T = " + T1 + "s (5/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]);

subplot(3, 1, 2);
plot(f, db2); grid on;
title("T = " + T1 + "s (15/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]);

subplot(3, 1, 3);
plot(f, db3); grid on;
title("T = " + T1 + "s (35/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]);

%% Part 4
% First, for my time step I selected an f0. This was an arbitrary selection
% but I assumed 100 Hz may make things easy since it is a round number. I
% set a value for the number of points N. I increased this until all of the
% plots had very smooth curves when zooming in to be sure that the
% correct features were captured. The range for t values was selected to be
% centered around 0 and give exactly 50 cycles of the sinusoid defined by
% f0 containing N points.

% The T values I selected to plot were based on the range given in the lab.
% I chose 5, 15, and 35 so there were values on each end and the middle of
% the acceptable range, and kept them as multiples of 5 for simplicity.
% Upon plotting, it was clear that these 3 would display trends well as
% there are clear differences between each of the plots.

% Finally, the frequency range was selected to be centered around f0. I
% chose lower and upper bounds that showed the peak well and some of the
% information around it. It did not seem necessary to show a very large
% range, since the plot appears to slowly asymptote, so I limited it to
% focus on the center frequency and some frequencies near it. I increased
% the resolution of f until all of the plots were very smooth. My original
% selection had some jagged edges on the 3rd T value, so I increased it
% until that plot was completely smooth and the shape looked as expected.

%% Part 5
% find the 1/2 of the max for the spectra U(f)
hmax1 = max(db1) - 6;
hmax2 = max(db2) - 6;
hmax3 = max(db3) - 6;

% use bool to find what range the max is achieved across
idx1 = (db1 >= hmax1);
idx2 = (db2 >= hmax2);
idx3 = (db3 >= hmax3);

% now get the frequency at the start and end
arr1 = f(idx1);
fl1 = arr1(1); fh1 = arr1(end);
arr2 = f(idx2);
fl2 = arr2(1); fh2 = arr2(end);
arr3 = f(idx3);
fl3 = arr3(1); fh3 = arr3(end);

% now calculate bandwidth
BW1 = fh1 - fl1;
BW2 = fh2 - fl2;
BW3 = fh3 - fl3;

% we can get the number of each cycles from T and determine the ratio
% between T and the 6dB bandwidth (after some experimentation)
r1 = BW1/f0*T1;
r2 = BW2/f0*T2;
r3 = BW3/f0*T3;

% we can see from the ratios that the fractional bandwidth is directly
% related to the number of pules. this creates a constant 0.012
% we can determine from this that the fractional bandwidth is related to
% the constant divided by the number of samples: BW/f0 = 0.012/T
% this result makes sense as we saw for a gaussian pulse, the bandwidth was
% related to a constant divided by the pulse width sigma, and this
% relationship is very similar.

% This relationship relates to how the Dirac delta and its Fourier
% transform function. A dirac delta is an incredibly strong pulse at
% only one location, which gives a Fourier transform of a horizontal line
% at 1. As you spread out from a delta function, the fourier transform will
% converge into a peak until you have created the opposite (one frequency
% peak, horizontal line in time domain). This is essentially the
% relationship we have determined. As the number of pulses goes up, the
% bandwidth will get smaller and smaller. As this continues, the bandwidth
% will approach 0 as T goes to infinity. This is just like the delta
% behaviour, except instead of having a perfect horizontal we have an
% oscillating pulse. The opposite happens as T goes to 0. The bandwidth will
% approach infinity. The relationship we determined between the bandwidth
% and the number of pulses follows the same concept as the dirac delta and
% its fourier transform.

% for fun, plot the bandwidths
% create figure
figure;

% plot each subplot
subplot(3, 1, 1);
hold on; plot(f, db1); xline([fl1 fh1], color="red"); hold off; grid on;
title("T = " + T1 + "s (5/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]); legend("Signal", "6dB BW");
subplot(3, 1, 2);
hold on; plot(f, db2); xline([fl2 fh2], color="red"); hold off; grid on;
title("T = " + T2 + "s (15/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]); legend("Signal", "6dB BW");
subplot(3, 1, 3);
hold on; plot(f, db3); xline([fl3 fh3], color="red"); hold off; grid on;
title("T = " + T3 + "s (35/f0)"); ylabel("Amplitude (dB)"); xlabel("Frequency (Hz)");
ylim([0 100]); xlim([f(1) f(end)]); legend("Signal", "6dB BW");