%%
% Slade Brooks
% brooksl@mail.uc.edu
% 11.05.25
% BME6013C
% Lab 09

clear variables
close all

%% Part 1
% read in original image
im = imread("p_6_1.jpg");

% convert to double
im = double(im);

% read in noisy image
data = load("p_3_4.mat");
im2 = data.noisy_MRI_image;

% set up range based on size (step of 0.1 cm)
dx = 0.1; dy = dx;
N = size(im, 1);
x = -N*dx:dx:N*dx-dx;
y = x;

% display images
figure
% plot original image
subplot(1, 2, 1);
imagesc(x, y, im);
colormap("gray"); axis image;
title("Original Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

% plot noisy image
subplot(1, 2, 2);
imagesc(x, y, im2);
colormap("gray"); axis image;
title("Noisy Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

%% Part 2
% edit example for setting up DFT
% Define a spatial-frequency grid as well, with zero
% spatial frequency at center (for plot after fftshift)
umax = 1/(2*dx); % Maximum spatial frequency for DFT, cm^-1
du = 1/(N*dx); % Spatial frequency step for DFT, cm^-1

% Corresponding vectors of spatial-frequency values
u = -umax:du:umax-du;
v = u;

% do DFT (and shift) of both images
IM = fftshift(abs(fft2(im)));
IM2 = fftshift(abs(fft2(im2)));

% normalize
IM = IM/max(IM, [], "all");
IM2 = IM2/max(IM2, [], "all");

% set limits for plotting
dBrange = 50; % Decibel range for DFT display

% plot both images
figure

% original
subplot(1, 2, 1);
imagesc(u, v, 20*log10(IM));
colormap("gray"); axis image;
title("DFT of Original Image"); xlabel("u (1/cm)"); ylabel("v (1/cm)");
xlim([u(1) u(end)]); ylim([v(1) v(end)]); c = colorbar; c.Label.String = "Amplitude (dB)";
clim([-dBrange 0]);

% noisy
subplot(1, 2, 2);
imagesc(u, v, 20*log10(IM2));
colormap("gray"); axis image;
title("DFT of Noisy Image"); xlabel("u (1/cm)"); ylabel("v (1/cm)");
xlim([u(1) u(end)]); ylim([v(1) v(end)]); c = colorbar; c.Label.String = "Amplitude (dB)";
clim([-dBrange 0]);

%% Part 3
% create filter following example
[fx, fy] = meshgrid(u, v);
fmag_2D = sqrt(fx.^2 + fy.^2);
alpha = 0.07/dx; % Gaussian width, normalized to sampling rate
filter_lp_2D = fftshift(exp(-fmag_2D.^2/(2*alpha^2)));
im2_lp = ifft2(fft2(im2).*filter_lp_2D);

% plot original and filtered
figure

% original
subplot(1, 2, 1);
imagesc(x, y, im2);
colormap("gray"); axis image;
title("Noisy Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

% filtered
subplot(1, 2, 2);
imagesc(x, y, im2_lp);
colormap("gray"); axis image;
title("Filtered Noisy Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

%% Part 4
% set up for loop
alphas = 0:0.005:0.5/dx;
err = zeros(size(alphas));
i = 1;

% loop through filter params
for a = 0:0.005:0.5/dx
    % calculate filter like example in part 3
    alpha = a; % Gaussian width, normalized to sampling rate
    filter_lp_2D = fftshift(exp(-fmag_2D.^2/(2*alpha^2)));

    % apply filter to image
    im2_filt = ifft2(fft2(im2).*filter_lp_2D);

    % get difference btn original and filtered
    diff = im2_filt - im;

    % calculate error ratio
    err(i) = rms(diff, "all")/rms(im, "all");
    i = i + 1;
end

% plot results
figure
hold on
plot(alphas*dx, err)
title("Error Ratio for Filter Strength")
xlabel("\alpha/f_s"); xlim([alphas(1)*dx alphas(end)*dx]);
ylabel("error ratio"); ylim([0 1]);
grid on

% grab minimum
[errmin, Lmini] = min(err);
Lmin = alphas(Lmini);
disp("minimum error is " + errmin + " at alpha=" + Lmin*dx)
scatter(Lmin*dx, errmin); hold off

% display final images
figure
% plot original image
subplot(1, 2, 1);
imagesc(x, y, im);
colormap("gray"); axis image;
title("Original Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

% plot filtered noisy image
% calculate filter like example in part 3
alpha = Lmin; % Gaussian width, normalized to sampling rate
filter_lp_2D = fftshift(exp(-fmag_2D.^2/(2*alpha^2)));

% apply filter to image
im2_filt = ifft2(fft2(im2).*filter_lp_2D);
subplot(1, 2, 2);
imagesc(x, y, im2_filt);
colormap("gray"); axis image;
title("Filtered Noisy Image"); xlabel("x (cm)"); ylabel("y (cm)");
xlim([x(1) x(end)]); ylim([y(1) y(end)]); c = colorbar; c.Label.String = "Amplitude";
clim([min(min(im, [], "all"), min(im2, [], "all")) max(max(im, [], "all"), max(im2, [], "all"))]);

% display final image FFTs
IM2_filt = fftshift(abs(fft2(im2_filt)));

% normalize
IM2_filt = IM2_filt/max(IM2_filt, [], "all");

% set limits for plotting
dBrange = 50; % Decibel range for DFT display

% plot both images
figure

% original
subplot(1, 2, 1);
imagesc(u, v, 20*log10(IM));
colormap("gray"); axis image;
title("DFT of Original Image"); xlabel("u (1/cm)"); ylabel("v (1/cm)");
xlim([u(1) u(end)]); ylim([v(1) v(end)]); c = colorbar; c.Label.String = "Amplitude (dB)";
clim([-dBrange 0]);

% noisy
subplot(1, 2, 2);
imagesc(u, v, 20*log10(IM2_filt));
colormap("gray"); axis image;
title("DFT of Filtered Noisy Image"); xlabel("u (1/cm)"); ylabel("v (1/cm)");
xlim([u(1) u(end)]); ylim([v(1) v(end)]); c = colorbar; c.Label.String = "Amplitude (dB)";
clim([-dBrange 0]);

%% Part 5
% The minimum error was at a value of alpha/fs=0.198. This means for our dx
% of 0.1 cm, the actual value of alpha is 1.98 1/cm. This is very
% interesting when looking at the DFT of the images. In the original image,
% it appears that there is no significant frequency content outside of |u|,
% |v| = 2 1/cm. In the DFT of the noisy image, there is frequency content
% visible everywhere. It makes a lot of sense then that the most accurate
% filter removes content outside of 1.98 1/cm to better match the original
% image. We can see from the final DFT of the original and the best filter
% that the original has frequency content in a square from roughly -2 to 2
% u and v, whereas the DFT of the best filter has a circle from
% roughly -2 to 2 u and v. The filtering does not perfectly recreate the
% frequency content of the original, but it does appear to significantly
% reduce the higher frequency noise that was visible initially.