clearvars; clc;
close all;

T = readtable('lab3_002_ALL');
data = table2array(T);

% h = figure();
% set(h,'WindowStyle','docked');

time = data(6:end, 1);
wave = data(6:end, 2);



% x = eyediagram(wave(7:end), 15, 6);

h = figure();
set(h,'WindowStyle','docked');
plot(time, wave);



ss = 40e-12; %Sample Period
fs = 1/ss; %Sample Frequency
L = 2^12; %Sampling Length
tt = (-L/2:L/2-1) .* ss; % Time Vector

wave_p = wave(1:L)

% 200 ns
% 2 GHz = .5 ns  - 500ps, 400 pulses in 200ns
% 5e3pts, 40ps/pt, 12.5 pts per pulse,  
% interp - 10e3 pts, 20ps/pt, 25pts per pulse



% Define Frequency Points
ff = (-L/2:L/2-1)/(ss*L);
ghz = 1e9;
fghz = ff/ghz;
% 
% % Numerical FFT
% fft_E = fftshift(fft(E));
% fft_E_norm = abs(fft_E).*ss;
% fft_E_unity = fft_E_norm/(tau*sqrt(pi));

psd = (abs(fftshift(fft(wave_p))).*ss).^2

plot(fghz, psd); hold on;
xlim([-5, 5]);
xticks(linspace(-5, 5, 11));
xlabel("Frequency [GHz]");
ylabel("PSD");
title("PSD of 2 GHz Signal and 4 GHz Scope Bandwidth");

print -depsc psd