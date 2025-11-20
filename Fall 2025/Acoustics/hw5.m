clear variables; close all;

% define properties from appendix
rhoSW = 1026;
cSW = 1500;
ZSW = 1.54e6;
rhoQS = 2070;
cQS = 1730;
ZQS = 3.58e6;

%% Part A
% set range of angles to check
thetai = -90:1:90;

% calculate xmission angle for each
thetat = asind(cQS/cSW.*sind(thetai));

% plot results
figure
plot(thetai, thetat)
xlim([thetai(1) thetai(end)]); xlabel("Incident Angle (deg)")
ylim([-100 100]); ylabel("Transmission Angle (deg)")
grid on

%% Part B
% calc critical angle
thetac = asind(cSW/cQS);

%% Part C
% calc reflection coef (magnitude)
Rbar = abs((ZQS./cosd(thetat) - ZSW./cosd(thetai))./(ZQS./cosd(thetat) + ZSW./cosd(thetai)));

% calc xmission coef (magnitude)
Tbar = abs(2*(ZQS./cosd(thetat))./(ZQS./cosd(thetat) + ZSW./cosd(thetai)));

% plot results
figure
hold on
plot(thetai, Rbar); plot(thetai, Tbar)
hold off
xlim([thetai(1) thetai(end)]); xlabel("Incident Angle (deg)")
ylim([0 2.5]); ylabel("Pressure Coefficient")
legend(["Reflected", "Transmitted"])
grid on

%% Part D
% calculate reflection intensity coeff
RI = Rbar.^2;

% plot
figure
plot(thetai, RI)
xlim([thetai(1) thetai(end)]); xlabel("Incident Angle (deg)")
ylim([0 1.5]); ylabel("Reflected Intensity Coefficient")
grid on