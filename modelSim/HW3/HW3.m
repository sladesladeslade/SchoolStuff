% Slade Brooks
% Model/Sim HW3

% make transfer fxns from hand calc ones
a = tf(15, [0, 5, 7]);
b = tf(5, [3, 30, 63]);
c = tf(4, [1, 4, 21]);

% calc step fxn for each
[Aa, ta] = step(a);
[Ab, tb] = step(b);
[Ac, tc] = step(c);

% plot step fxns
figure(1);
plot(ta, Aa);
ylabel("Amplitude")
xlabel("Time (s)")
title("Step Response (a)")
figure(2);
plot(tb, Ab);
ylabel("Amplitude")
xlabel("Time (s)")
title("Step Response (b)")
figure(3);
plot(tc, Ac);
ylabel("Amplitude")
xlabel("Time (s)")
title("Step Response (c)")