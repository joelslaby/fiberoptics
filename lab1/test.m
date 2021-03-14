close all; clearvars; clc;

Fs = 1e11;
delx = 1/Fs;

tb = 300e-9;
t = -tb:1/Fs:tb;
N = length(t);

tau_fwhm = 100e-12;
tau = tau_fwhm/sqrt(2*log(2));
Eo = 1;

I = Eo^2 .* exp(-2.*t.^2./tau^2);
plot(t, abs(I))
xlim([-300e-12, 300e-12])

phi = fft(I)/N;
phi = fftshift(phi);
delk = 2*pi/(N*delx);
k = (-N/2:N/2-1)*delk;
figure();
plot(k, abs(phi))