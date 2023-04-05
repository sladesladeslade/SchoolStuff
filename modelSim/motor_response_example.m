%% Written by Prof. Donghoon Kim

clear all; close all; clc

 

%% Example 6.6.2

KT = 0.05;  % Motor constant (N-m/A)

Kb = KT;    % Motor constant (N-m/A)

La = 2e-3;  % Inductor (H)

Ra = 0.5;   % Resistor (Ohm)

I = 9e-5;   % Inertia (kg-m^2)

c = 1e-4;   % (N-m-s/rad)

 

% current

current_tf = tf([I,c], [La*I, Ra*I+c*La, c*Ra+Kb*KT]);

 

% speed:

speed_tf = tf(KT, [La*I, Ra*I+c*La, c*Ra+Kb*KT]);

 

% Step voltage

[current, tc] = step( current_tf );

[speed, ts] = step( speed_tf );

figure(1);

subplot(2,1,1), plot(tc,10*current), grid, hold, xlabel('t (s)'), ylabel('Current (A)')

subplot(2,1,2), plot(ts,10*speed), grid, hold, xlabel('t (s)'), ylabel('Speed (rad/s)')

 

% Modified step input

t = (0:0.0001:0.07);

v = 10*(1 - exp( -t/0.01 ));

% lsim(SYS,U,T) plots the time response of the dynamic system SYS to the

% input signal described by U and T.

ia = lsim( current_tf, v, t );

omega = lsim( speed_tf, v, t );

%figure; lsim( current_tf, v, t );

figure, plot(t, ia, t, v), grid;

xlabel('t (s)'), ylabel('Amplitude'), legend('i_a', 'v')

figure(1);

subplot(2,1,1), plot(t,ia), xlim([t(1) t(end)])

subplot(2,1,2), plot(t,omega), xlim([t(1) t(end)])