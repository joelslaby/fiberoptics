close all; clearvars; clc;

%Gaussian pulse centered at t=0 with 100ps FWHM

tau_fwhm = 100e-12;
tau = tau_fwhm/sqrt(2*log(2));
Eo = 1;

tb = 300e-10;
Fs = 1e12; 
t = -tb:1/Fs:tb;
L = length(t);
n = 2^nextpow2(L);




E = Eo .* exp(-t.^2./tau^2);
I = abs(E.^2);
dbI = 10*log10(I);

figure();
t1 = t/1e-12;
plot(t1, I, 'Linewidth', 2);
xlim([-300, 300]);
xlabel("Time [ps]");
ylabel("Magnitude |E(t)|^2");
title("Gaussian Pulse");
print -depsc timeGauss

figure();
% plot(t1, dbI);
semilogy(t1, I, 'Linewidth', 2);
ylim([1e-3, 1]);
xlabel("Time [ps]");
ylabel("Magnitude |E(t)|^2");
title("Semi-log Time Plot of Gaussian Pulse")
print -depsc semilogTime

Y = fftshift(fft(E, n));


% my = abs(Y.^2)/n;
f = (-n/2:n/2-1)*Fs/n;

figure();
f1 = f/1e9;
psd_g = abs(Y/n).^2
norm_coef = (psd_g(find(f==0)));
psd_norm = psd_g./norm_coef;
plot(f1, psd_norm, 'Linewidth', 2);
xlim([-10, 10]);
title("Frequency Response of Gaussian Pulse")
xlabel("Frequency [GHz]");
ylabel("PSD");
print -depsc freqGauss


% [pxx,w] = 
% plot(w,10*log10(pxx))
% [pxx, w] = periodogram(I,rectwin(length(I)),length(I),Fs);
% w1 = w/1e12;
% periodogram(I,rectwin(length(I)),length(I),Fs);
% plot(w1,10*log10(pxx))
% f = Fs*[(0:(n))/n];
% P = abs(Y/n).^2;
% P = [P((n/2):(n)), P(1:(n/2))]
% 
% plot(f,P) 
% title('Gaussian Pulse in Frequency Domain')
% xlabel('Frequency (f)')
% ylabel('|P(f)|^2')


