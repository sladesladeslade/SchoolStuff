%%
% Slade Brooks
% brooksl@mail.uc.edu
% 10.15.25
% BME6013C
% Lab 07

clear variables
close all

%% Part 1
% steal spatial grid from example code (switch to cm)
% Define a spatial grid for all examples
xmax = 4; dx = 0.01; % spatial range and step (cm)
ymax = xmax; dy = dx; % vertical dimensions matching horizontal
xvec = -xmax:dx:xmax;
yvec = -ymax:dy:ymax;
[x,y] = meshgrid(xvec,yvec);

% Define a spatial-frequency grid as well (switch to 1cycle/cm max)
umax = 1; du = 0.004;
vmax = umax; dv = du;
uvec = -umax:du:umax;
vvec = -vmax:dv:vmax;
[u,v] = meshgrid(uvec,vvec);

% steal signal and fourier xform from example, update sigma to be both
% directions
sigmax = 1.5; sigmay = sigmax;
g = exp(-x.^2/(2*sigmax^2)).*exp(-y.^2/(2*sigmay^2))/(2*pi*sigmax*sigmay);
G = exp(-2*pi^2*(u.^2*sigmax^2 + v.^2*sigmay^2));

% plot signal
figure
imagesc(xvec, yvec, g, [-max(g,[],"all") max(g,[],"all")])
c = colorbar; c.Label.String = "Amplitude (cm)";
xlabel("x (cm)"); xlim([xvec(1) xvec(end)]);
ylabel("y (cm)"); ylim([yvec(1) yvec(end)]);
title("g(x, y)")
axis square

% plot F.T.
figure
imagesc(uvec, vvec, G, [-1 1])
c = colorbar; c.Label.String = "Amplitude (cm)";
xlabel("u (cycle/cm)"); xlim([uvec(1) uvec(end)]);
ylabel("v (cycle/cm)"); ylim([vvec(1) vvec(end)]);
title("G(u, v)")
axis square

%% Part 2
% create array of L values to try (cm)
L = 0:0.01:10;

% now loop through L values
err = zeros(size(L));
for i = 1:length(L)
    % grab current L to test
    Lx = L(i); Ly = Lx;

    % steal F.T. of box fxn
    B = sinc(u*Lx).*sinc(v*Ly);

    % get difference btn box and gauss
    diff = B - G;

    % calculate error ratio
    err(i) = rms(B - G, "all")/rms(G, "all");
end

% plot results
figure
hold on
plot(L, err)
title("error ratio for box function F.T.")
xlabel("L (cm)"); xlim([L(1) L(end)]);
ylabel("error ratio"); ylim([0 max(err)]);
grid on

% grab minimum
[errmin, Lmini] = min(err);
Lmin = L(Lmini);
disp("minimum error is " + errmin + " at L=" + Lmin + " cm")
scatter(Lmin, errmin); hold off

%% Part 3
% use the found L to minimize error
Lopt = Lmin;

% steal signal and f.t. from example
Lx = Lopt; Ly = Lx;
b = ((abs(x)<=Lx/2)&(abs(y)<=Ly/2)) / (Lx*Ly);
B = sinc(u*Lx).*sinc(v*Ly);

% plot signal
figure
imagesc(xvec, yvec, b, [-max(g,[],"all") max(g,[],"all")])
c = colorbar; c.Label.String = "Amplitude (cm)";
xlabel("x (cm)"); xlim([xvec(1) xvec(end)]);
ylabel("y (cm)"); ylim([yvec(1) yvec(end)]);
title("b(x, y)")
axis square

% plot F.T.
figure
imagesc(uvec, vvec, B, [-1 1])
c = colorbar; c.Label.String = "Amplitude (cm)";
xlabel("u (cycle/cm)"); xlim([uvec(1) uvec(end)]);
ylabel("v (cycle/cm)"); ylim([vvec(1) vvec(end)]);
title("B(u, v)")
axis square

%% Part 4
% plot fxn x-secs
figure; hold on;
plot(xvec, g(:,yvec==0)); plot(xvec, b(:,yvec==0)); hold off;
xlabel("x (cm)"); xlim([xvec(1) xvec(end)]);
ylabel("Amplitude (cm)"); ylim([-max(g,[],"all") max(g,[],"all")]);
title("g(x, 0) and b(x, 0)")
grid on
legend(["g" "b"])

% plot F.T. x-secs
figure; hold on;
plot(uvec, G(:,vvec==0)); plot(uvec, B(:,vvec==0)); hold off;
xlabel("u (cycles/cm)"); xlim([uvec(1) uvec(end)]);
ylabel("Amplitude (cm)"); ylim([-1 1]);
title("G(u, 0) and B(u, 0)")
grid on
legend(["G" "B"])

%% Part 5
% We can see in the spatial domain, g and b are not very similar. They are
% both symmetrical in x and y. It appears the square would roughly cover
% the part of the plot of the gaussian that is not approaching zero, but
% the obvious difference is the square has all one value slightly below the
% maximum amplitude of the gaussian and does not have the circular shape.
% This difference is also clear in Figure 6. We can see in a slice through
% the center that the square is 0 except across the length, and reaches a
% little more than half of the magnitude of the gaussian. The gaussian is a
% smooth curve that approaches 0 on the ends.

% In the frequency domain, the plots are much more similar. The magnitude,
% shape, and size of the narrow pulse at 0, 0 are almost identical. The
% F.T. of the box function is a very good approximation of the F.T. of the
% gaussian pulse. However, there are some oscillations visible in the box
% function that are not present in the gaussian. The gaussian has a
% magnitude of 0 everywhere outside of the main pulse, but the box function
% oscillates around 0. It also oscillates mostly in the x and y directions
% only, with only some small oscillation at different angles, so there is
% no radial symmetry. This is very visible in the slice in Figure 7, where
% we can see the box does a good job approximating the central pulse and
% its magnitude, but oscillates around 0 instead of being 0 outside the
% pulse. The oscillations appear to be weakening, so the function is likely
% much closer to 0 very far away, but there are significant oscillations
% near the pulse.

% I expect an image blurred by the gaussian would be very uniformly blurred
% and be very smooth. I think an image blurred with the box function would
% have a similar look, but I would expect some artifacts from the
% oscillation. I think the image would still be very smooth, but I expect
% the oscillations would appear and there would be waviness horziontally
% and vertically across the image. It would likely have the same smoothness
% but have some artifacting from the oscillation causing it to be distorted
% more than the gaussian would and possibly even look stripy.