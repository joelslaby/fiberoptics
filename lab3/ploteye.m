clearvars; clc; close all;

T = readtable('\media\lab3_000_ALL');
data = table2array(T);

h = figure();
set(h,'WindowStyle','docked');

time = data(7:end, 1);
wave = data(7:end, 2);

plot(time, wave);

% x = eyediagram(wave(7:end), 15, 6);

h = figure();
set(h,'WindowStyle','docked');

time_interp = linspace(time(1), time(end), length(time)*2-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

% time_eye = linspace();
offset = 20;
packet = 25;
for i=0:900
    plot(wave_interp((i*packet+offset):(((i+1)*packet+offset)-1))); hold on;
    ylim([-.5, .5])
end

% 100 ns
% 2 GHz = .5 ns  - 500ps, 8
% 25e3 40ps/pt