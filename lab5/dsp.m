clearvars; clc;
close all;

% SETUP
% 0 = No plots, Calculate param 
% 1 = Plot Gradient Eye Only, No param calculations, Used to center graph
% 2 = Plot everything, Calculate all param
% 3 = Plot all traces in Grid
plot_data = 2;

% 0-no fiber, no mod
% 1-no fiber, mod
% 2-fiber, mod
% 3-lab4 day1, wrong vPI
wave_type = 0;
select = 1;

filename = 'lab4/media_day2/lab4_009_ALL.csv';
    
T = readtable(filename);
data = table2array(T);

time = data(6:end, 1);
wave = data(6:end, 3);

ps = 1e-12;
mV = 1e-3;

h = figure();
set(h,'WindowStyle','docked');
% subplot(2, 1, 1);
endNbr = 400;
plot(time(1:endNbr)/ps, wave(1:endNbr)/mV); hold on;
xlim([time(1)/ps, time(endNbr)/ps]);

diff = [-1, 1];

waveDiff = conv(wave, diff);

% subplot(2, 1, 2);
% plot(time(1:endNbr)/ps, waveDiff(1:endNbr)/mV);

% pk = max(waveDiff);
% waveDiff(waveDiff>.7*pk) = 1;
% waveBin = waveDiff > .6*max(waveDiff) =
plot(time(1:endNbr)/ps, waveDiff(1:endNbr)/mV);
% findpeaks(waveDiff(1:endNbr)/mV);

% plot(time, )

