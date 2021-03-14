close all; clearvars; clc;

%Gaussian pulse centered at t=0 with 100ps FWHM
baud = 32e9;
tau_fwhm = 1/baud;
tau = tau_fwhm/sqrt(2*log(2));
Eo = 1;

tb = 300e-10;
Fs = 1e12; 
t = -tb:1/Fs:tb;
L = length(t);
n = 2^nextpow2(L);

E = Eo .* exp(-t.^2./tau^2);
I = abs(E).^2;
dbI = 10*log10(I);

% subplot(2, 2, 1);
t1 = t/1e-12;
% plot(t1, I);
% xlim([-100, 100]);
% xlabel("Time [ps]");
% ylabel("Magnitude |E(t)|^2");
% title("Gaussian Pulse")
% 
% subplot(2, 2, 2);
% plot(t1, dbI);
% semilogy(t1, I);
% ylim([1e-3, 1]);
% xlabel("Time [ps]");
% ylabel("Magnitude |E(t)|^2");
% title("Semi-log Time Plot of Gaussian Pulse");

Y = fftshift(fft(I, n));
lambda = 1550; %nm
D = 17; %ps/nm-km
c = 3e8;
unit_factor = 1e-21;
beta2 = D*lambda^2/(2*pi*c)*unit_factor;

f = (-n/2:n/2-1)*Fs/n;
z = [10, 100, 1000];
for i = 1:3
    prop(i, :) = exp(j*beta2*z(i)*(2*pi*f).^2/2);
end

f1 = f/1e9;
psd_g = abs(Y/n).^2;
norm_coef = (psd_g(find(f==0)));
psd_norm = psd_g./norm_coef;

for i = 1:3
    disperse(i, :) = prop(i, :).*Y;
    psd_disp(i, :) = abs(disperse(i, :)/n).^2;
%     norm_coef = (psd_disp(i,find(f==0)));
%     norm_coef = 1;
    psd_d_norm(i, :) = psd_disp(i, :)./norm_coef;
end

figure();
plot(f1, psd_norm, f1, psd_d_norm(1, :), '--','linewidth', 2);
xlim([-30, 30]);
title("Frequency Response of Gaussian Pulse")
xlabel("Frequency [GHz]");
ylabel("PSD");
legend("Initial Pulse", "Dispersed Pulse")
% print -depsc freqDisperse


% To = tau_fwhm/(2*sqrt(log(2)));
% ff = exp(-(2*pi*f).^2*To^2/2);
% plot(f1, psd_g);
% xlim([-30, 30]);
% title("Frequency Response of Gaussian Pulse");
% xlabel("Frequency [GHz]");
% ylabel("PSD");
% legend("psd_g", "psd_norm", "ideal");

% figure();
% f1 = f/1e9;
% plot(f1, angle(Y), f1, angle(disperse))
% xlim([-30, 30]);

figure();
xx = ifft(fft(E.^2, n));

for i=1:3
    yy(i, :) = ifft(disperse(i, :), n);
end

t = t/1e-12;
% subplot(1, 2, 1);
% plot(t, I, t, abs(xx(1:L)), '--')
% xlim([-100, 100]);
% legend("before", "after fft and ifft")
% subplot(1, 2, 2);
% plot(t, abs(xx(1:L)), t, I, t, abs(yy(1:L)));
hold on;
plot(t, abs(xx(1:L)), 'linewidth', 1.5);
lineVec = {'-.', '--', ':'};
for i=1:3
    plot(t, abs(yy(i, 1:L)), 'LineStyle', lineVec{i}, 'linewidth', 1.5);
end
xlim([-100, 100]);
legend("original", "10km", "100km", "1000km");
title("Gausian Pulse Dispersion over Varying Lengths");
xlabel("Time [ps]");
ylabel("Magnitude |E(t)|^2");
% print -depsc timeDisperse
