close all; clearvars; clc;

k = 1.380649e-23;
T = 300;
RL = 50;
delF = 12.5e9;
R = 0.59;
PdB = linspace(-20, 10, 100);
P = 10.^(PdB./10) * 1e-3;
NF = 3;
Fn = 10^(NF/10);

iT = sqrt(4 * k * T / RL *delF);

SNR = (R.*P).^2 ./ (4 * k * T / RL * delF * Fn);

SNRdB = 10.*log10(SNR);

h = figure();
set(h,'WindowStyle','docked');

plot(PdB, SNRdB, "linewidth", 2);

title("SNR vs. Average Received Power");
xlabel("Average Received Power [dBm]");
ylabel("SNR [dB]");

% print -depsc SNRvsRP

NEP = sqrt(4 * k * T * Fn / (RL * R^2));

% lambda = linspace(.9, 2.1, 25);
% R = [.25, .3, .34, .38, .41, .45, .48, .52, .55, .57, .58, .59, .59, .59, .57, .57, .6, ...
%     .68, .82, 1, 1.08, 1.05, .95, .65, .3];

% plot(lambda, R);

eta = 1.24*0.59/1.55;


RINdB = [-130:-5:-155];
RIN = 10.^(RINdB./10);

RINvar = zeros(length(RINdB), length(P));

for i=1:length(RINdB)
    RINvar(i, :) = P.^2 .* RIN(i) * delF;
end

h = figure();
set(h,'WindowStyle','docked');

semilogy(PdB, RINvar, "linewidth", 2);
legend("-130 dB/Hz", "-135 dB/Hz", "-140 dB/Hz", "-145 dB/Hz", "-150 dB/Hz", "-155 dB/Hz", 'Location', 'northwest');
title("RIN Variance vs. Received Optical Power");
xlabel("Average Received Optical Power [dBm]");
ylabel("RIN Variance [A^2]");

print -depsc RINvar