close all; clc; clearvars;

Q = linspace(0, 8, 80000);

ber = 1/2*erfc(Q/sqrt(2));
berapprox = exp(-Q.^2./2)./(Q.*sqrt(2*pi));
h = figure();
set(h,'WindowStyle','docked');

semilogy(Q, ber, Q, berapprox);
title("BER vs. Q");
xlabel("Q");
ylabel("BER");
yticks(logspace(-18, 3, 8));
ylim([1e-18, 1e3]);
legend("erfc function", "approximated erfc")

% print -depsc ber

Q(find(ber<1e-12,1))
Q(find(ber<1e-9,1))
Q(find(ber<1e-3,1))


b = 10e9; % Gbps
ber = [1e-9, 1e-12, 1e-151]
t = 10 * 1./ber * 1/b % each bit
% 1s | 16min 40s | 11 days 13 hours 46 min 40s

