%%
% Slade Brooks
% brooksl@mail.uc.edu
% 11.12.25
% BME6013C
% Lab 10

clear variables
close all

%% Part 1
% make random array
chan = 256;
N = 1024;
noise = randn(N, chan);

% display image
figure
imagesc(noise);
colormap("gray"); axis image;
title("Noise"); xlabel("Channel"); ylabel("t (s)");

%% Part 2
% make and apply window
Wi = rectwin(N);
noiseW = noise.*Wi;

% take DFT of each channel
Noise = fftshift(fft(noiseW, [], 1)/N);

% take average of magnitude^2 for each freq in different sets of channels
aveNoise14 = mean(abs(Noise(:,1:4)).^2, 2);
aveNoise116 = mean(abs(Noise(:,1:16)).^2, 2);
aveNoise164 = mean(abs(Noise(:,1:64)).^2, 2);
aveNoise1256 = mean(abs(Noise(:,1:256)).^2, 2);

% make shifted freq array for plotting
dx = 0.1;     % set a time step so we have frequencies
umax = 1/(2*dx); % Maximum spatial frequency for DFT, Hz
du = 1/(N*dx); % Spatial frequency step for DFT, Hz

% Corresponding vectors of spatial-frequency values
u = -umax:du:umax-du;

% plot each power spectra
figure
subplot(2, 2, 1)
plot(u, aveNoise14)
title("Estimated Power Spectra (Chan. 1-4)"); xlabel("Frequency (Hz)");
ylabel("Power (cm)"); grid on; ylim([0 3.5e-3]);
subplot(2, 2, 2)
plot(u, aveNoise116)
title("Estimated Power Spectra (Chan. 1-16)"); xlabel("Frequency (Hz)");
ylabel("Power (cm)"); grid on; ylim([0 3.5e-3]);
subplot(2, 2, 3)
plot(u, aveNoise164)
title("Estimated Power Spectra (Chan. 1-64)"); xlabel("Frequency (Hz)");
ylabel("Power (cm)"); grid on; ylim([0 3.5e-3]);
subplot(2, 2, 4)
plot(u, aveNoise1256)
title("Estimated Power Spectra (Chan. 1-256)"); xlabel("Frequency (Hz)");
ylabel("Power (cm)"); grid on; ylim([0 3.5e-3]);

% We can see that as we average more channels together, the power spectra
% seems to converge to a flat distribution. The averages of 4 and 16
% channels have a decent amount of peaks/valleys and a significant
% deviation from a mean of 1. As the number of samples increases, we can
% see the amount of deviation decrease and it clearly converges to a mean
% of 1. The plot of 256 channels is a small amount of noise around 1. It
% appears that the "structures" in the power spectra do remain as more
% channels are averaged, but the power in each structure is reduced which
% causes a clear convergence.

%% Part 3
% calculate standard deviation of each
sig14 = std(aveNoise14);
sig116 = std(aveNoise116);
sig164 = std(aveNoise164);
sig1256 = std(aveNoise1256);

% make a plot of it
figure
hold on

% add the line of 1/sqrt(Nave) scaled by N pts
Naves = 1:0.1:256;
plot(Naves, 1./sqrt(Naves)/N, "LineWidth", 3)

% plot each of the std calculated w/ the corresponding num of channels
scatter(4, sig14, 100, "filled"); scatter(16, sig116, 100, "filled");
scatter(64, sig164, 100, "filled"); scatter(256, sig1256, 100, "filled");
hold off
xlabel("\# of Chan. Averaged $(N_{ave}$)", Interpreter="latex"); ylabel("Standard Deviation");
title("Standard Deviation for # Channels Averaged"); grid on; legend(["$\frac{1/N}{\sqrt{N_{ave}}}$","4 Chan.","16 Chan.","64 Chan.","256 Chan."], Interpreter="latex", FontSize=12);
xlim([Naves(1) Naves(end)]); ylim([0 1e-3]);

%% Part 4
% do autocorr (multiply each fft by its own conjugate)
acN14 = fftshift(ifft(fftshift(aveNoise14.*conj(aveNoise14))));
acN116 = fftshift(ifft(fftshift(aveNoise116.*conj(aveNoise116))));
acN164 = fftshift(ifft(fftshift(aveNoise164.*conj(aveNoise164))));
acN1256 = fftshift(ifft(fftshift(aveNoise1256.*conj(aveNoise1256))));

% make time array (shifted)
ts = -N/2*dx:dx:N/2*dx - dx;

% plot results for each averaged signal
figure
subplot(2, 2, 1);
plot(ts, acN14);
title("IFFT{Autocor} (Chan. 1-4)"); xlabel("Lag (s)");
ylabel("Amplitude (cm)"); grid on; xlim([ts(1) ts(end)]); ylim([0 13e-7]);
subplot(2, 2, 2);
plot(ts, acN116);
title("IFFT{Autocor} (Chan. 1-16)"); xlabel("Lag (s)");
ylabel("Amplitude (cm)"); grid on; xlim([ts(1) ts(end)]); ylim([0 13e-7]);
subplot(2, 2, 3);
plot(ts, acN164);
title("IFFT{Autocor} (Chan. 1-64)"); xlabel("Lag (s)");
ylabel("Amplitude (cm)"); grid on; xlim([ts(1) ts(end)]); ylim([0 13e-7]);
subplot(2, 2, 4);
plot(ts, acN1256);
title("IFFT{Autocor} (Chan. 1-256)"); xlabel("Lag (s)");
ylabel("Amplitude (cm)"); grid on; xlim([ts(1) ts(end)]); ylim([0 13e-7]);

% The real autocorrelation of the gaussian white noise would be a single
% point at infinity at zero lag, with 0 everywhere else. We can see that
% all of the iffts are fairly similar to this. Obviously we can not expect
% a perfect delta recreation, but it is clear that as more samples are
% added, we get closer to a delta. The content everywhere around 0 should
% be 0, which it almost is for averaging all 256 channels. This effect
% makes sense as the F.T. of the signals got closer and closer to a flat
% line at 1, which would give an IFT of a perfect delta. The one main
% difference is that although the noise is reduced, the magnitude of the
% delta is also reduced for each increasing average. Ideally the magnitude
% of the zero component would increase closer to infinity as the noise is
% reduced, but we see the opposite effect.