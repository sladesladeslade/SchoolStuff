%%
% Slade Brooks
% brooksl@mail.uc.edu
% 11.19.25
% BME6013C
% Lab 11

clear variables
close all

%% read in both images
% set image url
imageURL = "https://media.istockphoto.com/id/465885010/photo/human-skull-x-ray-image.jpg?s=612x612&w=0&k=20&c=ePoyTJ-x1lgr5Gj2xG49HymniOX1ibOm4jX-CA20PPg=";

% read in image, convert to grayscale
I = im2double(max(webread(imageURL), [], 3));
AI = im2double(max(imread("jpeg.jpg"), [], 3));

% plot image to verify
figure();
subplot(1, 2, 1)
imshow(I); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 2, 2); imshow(AI); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;

% fix cropping of images to be more similar
I = I(10:500,25:end-25);
AI = AI(50:1059,:);
figure();
subplot(1, 2, 1)
imshow(I); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 2, 2); imshow(AI); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;

%% Part 3
%%% compute histograms of both images
figure()
hold on
histAI = histogram(AI, 256, "Normalization","probability");
histI = histogram(I, 256, "Normalization","probability");
hold off
legend(["AI Image" "Real Image"])
xlim([0 1]); grid on; %ylim([0 1]);
colormap(gray);
c = colorbar('southoutside'); c.TickLabels = []; c.Ticks = [];

% plot zoomed in version
figure()
hold on
histAI = histogram(AI, 256, "Normalization","probability");
histI = histogram(I, 256, "Normalization","probability");
hold off
legend(["AI Image" "Real Image"])
xlim([0 1]); grid on; ylim([0 0.05]);
colormap(gray);
c = colorbar('southoutside'); c.TickLabels = []; c.Ticks = [];

%%% compute spatial freq content of each image
% get size of images
NIx = size(I, 2);
NIy = size(I, 1);
NAIx = size(AI, 2);
NAIy = size(AI, 1);

% set physical dimension and correct
skW_cm = 18;       % Physical skull width in cm (example)
skW_px_AI = 914 - 16;     % Skull width in pixels for AI image
skW_px_I = 458 - 33;
pixelSize_cm_AI = skW_cm/skW_px_AI;
pixelSize_cm_I = skW_cm/skW_px_I;
dxAI = pixelSize_cm_AI;
dyAI = pixelSize_cm_AI;
dxI = pixelSize_cm_I;
dyI = pixelSize_cm_I;

% set up frequency vectors based on spatial
umaxAI = 1/(2*dxAI);
vmaxAI = 1/(2*dyAI);
umaxI = 1/(2*dxI);
vmaxI = 1/(2*dyI);
duAI = 1/(NAIx*dxAI);
dvAI = 1/(NAIy*dyAI);
duI = 1/(NIx*dxI);
dvI = 1/(NIy*dyI);
uAI = -umaxAI:duAI:umaxAI - duAI;
vAI = -vmaxAI:dvAI:vmaxAI - dvAI;
uI = -umaxI:duI:umaxI - duI;
vI = -vmaxI:dvI:vmaxI - dvI;

% now do the dft of a similar slice through the skull and the normalized 2D
% fft of each image
F_I = fftshift(abs(fft(I(125,:))))/max(abs(fft(I(125,:))), [], "all");
F_AI = fftshift(abs(fft(AI(245,:))))/max(abs(fft(AI(245,:))), [], "all");
F2_I = fftshift(abs(fft2(I)))/max(abs(fft2(I)), [], "all");
F2_AI = fftshift(abs(fft2(AI)))/max(abs(fft2(AI)), [], "all");

% plot FFTs
figure()
hold on
plot(uI*dxI, 20*log10(F_I));
plot(uAI*dxAI, 20*log10(F_AI)); grid on;
hold off
xlim([-0.5 0.5]); ylim([-80 0]);
xlabel("Spatial Frequency in x (u, cycles/cm)"); ylabel("Amplitude (dB)");
legend(["Real Image" "AI Image"]);

% plot 2D ffts
figure()
subplot(1, 2, 1); imagesc(uI*dxI, vI*dyI, 20*log10(F2_I)); axis image; colormap gray; clim([-80 0]); c = colorbar; c.Label.String = "Amplitude (dB)";
xlabel("Spatial Frequency in x * Pixel Size (u*dx)"); ylabel("Spatial Frequency in y * Pixel Size (v*dy)");
subplot(1, 2, 2); imagesc(uAI*dxAI, vAI*dyAI, 20*log10(F2_AI)); axis image; colormap gray; clim([-80 0]); c = colorbar; c.Label.String = "Amplitude (dB)";
xlabel("Spatial Frequency in x * Pixel Size (u*dx)"); ylabel("Spatial Frequency in y * Pixel Size (v*dy)");

%% Part 4
%%% first segment out tissue and remove
% tune to find tissue and set it as background
background = AI<20/255;
tissue = ~background & AI<125/255;
background = background | tissue;
edges = ~tissue & ~background & AI>225/255;
AI1 = AI;
AI1(background) = (0.997 + 1.994)/2/255;
AI1(edges) = 255/255;
rest = ~background & ~edges;

% spot correct weird artifact
AI1(909:1000,522:583) = 0;

% plot result
figure();
subplot(1, 3, 1)
imshow(I); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 3, 2)
imshow(AI); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 3, 3); imshow(AI1); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;

%%% gaussian low pass filter
% use matlab built-in 2d gaussian filter
sigma = 0.85;
AIfilt = imgaussfilt(AI1, sigma);

% plot result
figure();
subplot(1, 3, 1)
imshow(I); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 3, 2)
imshow(AI); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;
subplot(1, 3, 3); imshow(AIfilt); axis image;
xlabel("x (px)"); ylabel("y (px)"); colorbar;

% plot FFTs
figure()
hold on
F_AIfilt = fftshift(abs(fft(AIfilt(245,:))))/max(abs(fft(AIfilt(245,:))),[],"all");
plot(uI*dxI, 20.*log10(F_I));
plot(uAI*dxAI, 20.*log10(F_AIfilt));
plot(uAI*dxAI, 20.*log10(F_AI));
grid on;
hold off
xlim([-0.5 0.5]); ylim([-80 0]);
xlabel("Spatial Frequency in x (u, cycles/cm)"); ylabel("Amplitude (dB)");
legend(["Real Image" "Filtered AI Image" "Original AI Image"]);

% plot 2D ffts
F2_AIfilt = fftshift(abs(fft2(AIfilt)))/max(abs(fft2(AIfilt)),[],"all");
figure()
subplot(1, 3, 1); imagesc(uI*dxI, vI*dyI, 20*log10(F2_I)); axis image; colormap gray; clim([-80 0]); c = colorbar; c.Label.String = "Amplitude (dB)";
xlabel("Spatial Frequency in x * Pixel Size (u*dx)"); ylabel("Spatial Frequency in y * Pixel Size (v*dy)");
subplot(1, 3, 2); imagesc(uAI*dxAI, vAI*dyAI, 20*log10(F2_AI)); axis image; colormap gray; clim([-80 0]); c = colorbar; c.Label.String = "Amplitude (dB)";
xlabel("Spatial Frequency in x * Pixel Size (u*dx)"); ylabel("Spatial Frequency in y * Pixel Size (v*dy)");
subplot(1, 3, 3); imagesc(uAI*dxAI, vAI*dyAI, 20*log10(F2_AIfilt)); axis image; colormap gray; clim([-80 0]); c = colorbar; c.Label.String = "Amplitude (dB)";
xlabel("Spatial Frequency in x * Pixel Size (u*dx)"); ylabel("Spatial Frequency in y * Pixel Size (v*dy)");

% plot zoomed in version histogram
figure()
hold on
histAI = histogram(AI, 256, "Normalization","probability");
histI = histogram(I, 256, "Normalization","probability");
histAIfilt = histogram(AIfilt, 256, "Normalization","probability");
hold off
legend(["AI Image" "Real Image" "Filtered AI Image"])
xlim([0 1]); grid on; ylim([0 0.05]);
colormap(gray);
c = colorbar('southoutside'); c.TickLabels = []; c.Ticks = [];