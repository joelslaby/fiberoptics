clearvars; clc;
close all;

scale = 4;

T = readtable('lab3_006_ALL');
data = table2array(T);

% h = figure();
% set(h,'WindowStyle','docked');

time = data(6:end, 1);
wave = data(6:end, 2);

% plot(time, wave);

% x = eyediagram(wave(7:end), 15, 6);

h = figure();
set(h,'WindowStyle','docked');
subplot(3, 1, 1);

time_interp = linspace(time(1), time(end), length(time)*scale-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

offset = 13;
packet = 50;

time_eye = ((1:packet)-packet/2)*(20e-12/scale);

eye = reshape(wave_interp(offset:(packet*198+offset-1)), [packet, 198]);

ps = 1e-12;
time_eye_ps = time_eye/ps;
mV = 1e-3;
eye_mv = eye/mV;
plot(time_eye_ps, eye_mv);
% ylim([-.5, .5]);
xlim([-240*2/scale, 250*2/scale]);
title("Eye Diagram - 4 GHz signal - 4 GHz Bandwidth");
xlabel("Time [ps]");
ylabel("Voltage [V]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = readtable('lab3_007_ALL');
data = table2array(T);

% h = figure();
% set(h,'WindowStyle','docked');

time = data(6:end, 1);
wave = data(6:end, 2);

% plot(time, wave);

% x = eyediagram(wave(7:end), 15, 6);

subplot(3, 1, 2);

time_interp = linspace(time(1), time(end), length(time)*scale-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

time_eye = (1:25)*40e-12;

offset = 13;
packet = 50;

time_eye = ((1:packet)-packet/2)*(20e-12/scale);

eye = reshape(wave_interp(offset:(packet*198+offset-1)), [packet, 198]);

ps = 1e-12;
time_eye_ps = time_eye/ps;
mV = 1e-3;
eye_mv = eye/mV;
plot(time_eye_ps, eye_mv);
% ylim([-.5, .5]);
xlim([-240*2/scale, 250*2/scale]);
title("Eye Diagram - 4 GHz signal - 3 GHz Bandwidth");
xlabel("Time [ps]");
ylabel("Voltage [V]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = readtable('lab3_008_ALL');
data = table2array(T);

% h = figure();
% set(h,'WindowStyle','docked');

time = data(6:end, 1);
wave = data(6:end, 2);

% plot(time, wave);

% x = eyediagram(wave(7:end), 15, 6);

subplot(3, 1, 3);

time_interp = linspace(time(1), time(end), length(time)*scale-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

time_eye = (1:25)*40e-12;

offset = 13;
packet = 50;

time_eye = ((1:packet)-packet/2)*(20e-12/scale);

eye = reshape(wave_interp(offset:(packet*198+offset-1)), [packet, 198]);

ps = 1e-12;
time_eye_ps = time_eye/ps;
mV = 1e-3;
eye_mv = eye/mV;
plot(time_eye_ps, eye_mv);
% ylim([-.5, .5]);
xlim([-240*2/scale, 250*2/scale]);
title("Eye Diagram - 4 GHz signal - 2 GHz Bandwidth");
xlabel("Time [ps]");
ylabel("Voltage [V]");

print -depsc eyediag_compare_4gSig

