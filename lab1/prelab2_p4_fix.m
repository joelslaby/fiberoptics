clc; clearvars; close all;

h = figure();
set(h,'WindowStyle','docked');

baud = 32e9;
tau_fwhm = 1/baud;
tau = tau_fwhm/sqrt(2*log(2));
Eo = 1;
To = tau_fwhm/2*sqrt(log(2));

tb = 300e-10;
Fs = 1e12; 
t = -tb:1/Fs:tb;
L = length(t);
n = 2^nextpow2(L);

E = Eo .* exp(-t.^2./tau^2);
I = abs(E).^2;
dbI = 10*log10(I);

lambda = 1550; %nm
D = 17; %ps/nm-km
c = 3e8;
unit_factor = 1e-3;
beta2 = D*lambda^2/(2*pi*c)*unit_factor;

f = (-n/2:n/2-1)*Fs/n;
f_ghz = f/1e9;
z = [10, 100, 1000];
for i = 1:3
    prop(i, :) = exp(j*beta2*z(i)*(2*pi*f_ghz).^2/2);
end

t1 = t/1e-12;
plot(t1, E, t1, I);
xlim([-100, 100]);
xlabel("Time [ps]");
ylabel("Magnitude |E(t)|^2");
title("Gaussian Pulse");

X = fftshift(fft(E, n));

for i = 1:3
    disperse(i, :) = prop(i, :).*X;
end

Y = ifft(X, n);
for i=1:3
    Z(i, :) = ifft(disperse(i, :), n);
end

plot(t1, I, 'linewidth', 1.5); hold on;
lineVec = {':', '--', '-.'};
for i=1:3
   plot(t1, abs(Z(i, 1:L)).^2, 'LineStyle', lineVec{i}, 'linewidth', 1.5);
end

xlim([-100, 100]);
% ylim([0, .2])
xlabel("Time [ps]");
ylabel("Magnitude |E(t)|^2");
title("Gaussian Pulse with Dispersion Pulse Broadening");
legend("original", "10km", "100km", "1000km")

tauZ = sqrt(1+(beta2*z(3)/(To*1e12)^2)^2)*To*2*sqrt(log(2))

To*1e12
T1 = sqrt((To*1e12)^2+(beta2*1e3*z(1))^2)

print -depsc timeDisperse
