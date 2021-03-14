close all; clc; clearvars;


% 001 - SYS source
% 002 - 6km
% 005 - 25km

% 001 - SYS source
% 002 - 6km
% 005 - 25km


fileNbr = ["001", "002", "005"];

T = readtable('day2_000_001.csv');
data_src = table2array(T);

T = readtable('day2_000_002.csv');
data_6k = table2array(T);

h = figure();
set(h,'WindowStyle','docked');

time_src = data_src(:, 1);
wave_src = data_src(:, 2);

time_6k = data_6k(:, 1);
wave_6k = data_6k(:, 2);

% plot(time, wave);

start_src = find(time_src>-5e-9, 1);
final_src = find(time_src<5e-9, 1, "last");
clipTime_src = time_src(start_src:final_src);
clipWave_src = wave_src(start_src:final_src);

start_6k = find(time_6k>-5e-9, 1);
final_6k = find(time_6k<5e-9, 1, "last");
clipTime_6k = time_6k(start_6k:final_6k);
clipWave_6k = wave_6k(start_6k:final_6k);

clipTimeDense_src = linspace(clipTime_src(1), clipTime_src(end), length(clipTime_src)*100);
interpWave_src = interp1(clipTime_src, clipWave_src, clipTimeDense_src, 'spline');

clipTimeDense_6k = linspace(clipTime_6k(1), clipTime_6k(end), length(clipTime_6k)*100);
interpWave_6k = interp1(clipTime_6k, clipWave_6k, clipTimeDense_6k, 'spline');

h = figure();
set(h,'WindowStyle','docked');

ns = 1e-9;
mV = 1e-3;

% plot(clipTimeDense/ns, interpWave/mV, "Linewidth", 2);
% title("");
% xlabel("Time [ns]");
% ylabel("Voltage [V]");
% xlim([-3, 3]);

h = figure();
set(h,'WindowStyle','docked');

clipTimeDense_src = clipTimeDense_src - clipTimeDense_src(find(interpWave_src == max(interpWave_src)));
clipTimeDense_6k = clipTimeDense_6k - clipTimeDense_6k(find(interpWave_6k == max(interpWave_6k)));

plot(clipTimeDense_src/ns, interpWave_src/mV, "Linewidth", 2); hold on;
plot(clipTimeDense_6k/ns, interpWave_6k/mV, "Linewidth", 2);
title("");
xlabel("Time [ns]");
ylabel("Voltage [V]");
xlim([-.5, .5]);
