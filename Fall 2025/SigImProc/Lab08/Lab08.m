%%
% Slade Brooks
% brooksl@mail.uc.edu
% 10.22.25
% BME6013C
% Lab 08

clear variables
close all

%% Part 1
% read in data from file
T = readtable("p_10_1.xls");

% store one of the channels as g(t)
g = T.Chan4';

% create time vector
fs = 128;
Ts = 1/fs;
N = size(g, 2);
t = 0:Ts:N*Ts - Ts;

% plot the signal
figure
plot(t, g)
xlim([t(1) t(end)]); xlabel("Time (s)");
ylim([-15 20]); ylabel("Amplitude (cm)");
title("g(t)"); grid on

%% Part 2
% take fft of signal
G = fft(g);

% now shift the fft to be centered
Gshift = fftshift(G);

% create frequency array w/ 0 at N/2 + 1
fm = fs/N.*(-N/2:N/2 - 1);

% plot fft
figure
plot(fm, Gshift)
xlim([fm(1) fm(end)]); xlabel("Frequency (Hz)");
ylim([-300 1000]); ylabel("Amplitude");
title("G(f)"); grid on

%% Part 3
% create figure with subplots
figure

% need a new version of fm centered for the flipped guy
fm2 = fs/N.*(-N/2 + 1:N/2);

% plot the real of the conjugate and G flipped
subplot(2, 1, 1)
hold on
plot(fm, real(conj(Gshift)), linewidth=1.5); plot(fm2, real(flip(Gshift)), linestyle="--", linewidth=1.5);
hold off
xlim([fm(1) fm(end)]); xlabel("Frequency (Hz)");
ylim([-300 1000]); ylabel("Amplitude");
title("Re\{G(f)*\} and Re\{G(-f)\}"); grid on
legend(["Re\{G(f)*\}" "Re\{G(-f)\}"])

% plot the imaginary of the conjugate and the negative part of G
subplot(2, 1, 2)
hold on
plot(fm, imag(conj(Gshift)), linewidth=1.5); plot(fm2, imag(flip(Gshift)), linestyle="--", linewidth=1.5);
hold off
xlim([fm(1) fm(end)]); xlabel("Frequency (Hz)");
ylim([-500 500]); ylabel("Amplitude");
title("Im\{G(f)*\} and Im\{G(-f)\}"); grid on
legend(["Im\{G(f)*\}" "Im\{G(-f)\}"])

% These plots show symmetry because the signals overlap. G(-f) would be a
% flipped version of G(f). The real part of the conjugate of G(f) would
% just be the original signal. If a symmetrical signal was flipped, it
% would be identical to the original signal, so we are verifying that. For
% the imaginary part, taking the conjugate flips the imaginary part. Along
% the same lines, if this is the same as the imaginary part of the original
% signal flipped, then the signal must be symmetrical. We are proving
% through the conjugate and flipping the signal about the y axis that it is
% the same in both directions when flipped, meaning it must be symmetrical.

%% Part 4
% create rectangular window
Wi = rectwin(N);

% multiply time signal by rect window
grect = g.*Wi';

% create blackman window
Wii = blackman(N);

% multiply time signal by blackman window
gblack = g.*Wii';

% retake dft of both signals
Grect = fftshift(fft(grect));
Gblack = fftshift(fft(gblack));

% convert to dB
Grect = 20*log10(Grect);
Gblack = 20*log10(Gblack);

% plot signals
figure
hold on
plot(fm, Grect); plot(fm, Gblack)
hold off
xlim([fm(1) fm(end)]); xlabel("Frequency (Hz)");
ylim([-20 60]); ylabel("Amplitude");
title("G(f) with Rectangular and Blackman Window"); grid on
legend(["Rectangular Window" "Blackman Window"]);

%% Part 5
% It is clear that the blackman window removes some power from the original
% signal. The 0 frequency is reduced by almost 10 dB, and all of the power
% at frequencies above |20| Hz is reduced as well. The lower frequency
% power is similar in shape and magnitude to the rectangular, with some
% small differences. Any high frequency content is significantly reduced.
% The overall shape of the spectra stays the same, but there are some clear
% differences between the two. The blackman window tends to cause a much larger
% variation from the "mean" at the high frequencies, but also does remove a
% few very large peaks.

%% Part 6
% first set new frequency to get interpolation factor
fs2 = 1024;
scalefac = floor(fs2/fs);

% use factor to get new length of signal and arrays
N2 = N*scalefac;
f2 = (0:N2 - 1)/N2*fs2;
dt2 = Ts/scalefac; t2 = 0:dt2:N2*dt2 - dt2;

% create zero-padded fft and force symmetry
G2 = zeros(1, N2);
G2(1:N/2 + 1) = G(1:N/2+1);
G2(N2-N/2 + 2:N2) = G(N/2 + 2:N);
G2(N/2 + 1) = G2(N/2 + 1)/2;
G2(N2 - N/2 + 1) = G2(N/2 + 1);

% do inverse FFT on zero-padded fft
g2 = ifft(G2)*scalefac;

% plot old and new signal
figure
hold on
plot(t, g); plot(t2, g2);
hold off
xlim([0.5 1.0]); xlabel("Time (s)");
ylim([-15 20]); ylabel("Amplitude (cm)");
title("g(t)"); grid on; legend(["Original Signal" "Upsampled Signal"])

%% Part 7
% The interpolated signal is much smoother than the original. It appears to
% be identical to the original signal and clearly passes through all
% of the original data points. This would imply that the Nyquist frequency
% has been satisfied when sampling, otherwise I would expect to see
% artifacts of the upsampling in the new signal. If the signal wasn't
% sampled at a high enough rate, the zero-padding should cause incorrect
% upsampling that clearly does not match the original data. This also makes
% sense based on the power spectra. There is frequency content up to a
% magnitude of 63.5 Hz, meaning our Nyquist frequency would be 127 Hz. Our
% sample was taken at 128 Hz, so this should satisfy the requirement.