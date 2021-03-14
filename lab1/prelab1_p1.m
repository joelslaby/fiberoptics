clc; clearvars; close all;

%prelab 1 - p1
A = 7.81e11; %dB/km
air = 48.48; %um
B = .75; %dB-um^4/km

%a
lambda = linspace(.8, 1.65, 1e3);
alpha = A.*exp(-air./lambda)+B./lambda.^4;


semilogy(lambda, alpha); hold on;
xlim([.8, 1.65]);
ylim([.1, 10]);
title("Expected dB/km Fiber Loss over Desired Wavelengths")
xlabel("Wavelength [\mum]");
ylabel("Loss [dB/km]");
%b
lambda = [.850, 1.310, 1.550];
alpha = A.*exp(-air./lambda)+B./lambda.^4;
plot(lambda, alpha, '.black', 'MarkerSize', 15);

 print -depsc FiberLoss
