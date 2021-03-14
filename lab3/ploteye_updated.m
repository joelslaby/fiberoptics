clearvars; clc;
close all;

T = readtable('lab3_009_ALL');
data = table2array(T);

% h = figure();
% set(h,'WindowStyle','docked');

time = data(6:end, 1);
wave = data(6:end, 2);

% plot(time, wave);

% x = eyediagram(wave(7:end), 15, 6);

h = figure();
set(h,'WindowStyle','docked');


time_interp = linspace(time(1), time(end), length(time)*6-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

time_eye = (1:25)*40e-12;

offset = 13;
packet = 50;

time_eye = ((1:packet)-packet/2)*20e-12/6;

eye = reshape(wave_interp(offset:(packet*198+offset-1)), [packet, 198]);

% for i=0:198
%     subplot(2, 1, 1)
%     xx = wave_interp((i*packet+offset):(((i+1)*packet+offset)-1));
%     plot((50*i+1:50*i+50), wave_interp((i*packet+offset):(((i+1)*packet+offset)-1))); hold on;
%     ylim([-.5, .5])
%     
%     subplot(2, 1, 2)
%     plot(time_eye, wave_interp((i*packet+offset):(((i+1)*packet+offset)-1))); hold on;
    %plot(wave((i*packet+offset):(((i+1)*packet+offset)-1))); hold on;
    
% end

ps = 1e-12;
time_eye_ps = time_eye/ps;
mV = 1e-3;
eye_mv = eye/mV;
plot(time_eye_ps, eye_mv);
ylim([-400, 400]);
xlim([-240*2/2, 250*2/2]);
title("Eye Diagram - 6GHz signal - 4 GHz Bandwidth");
xlabel("Time [ps]");
ylabel("Voltage [V]");
print -depsc eyediag_6gSig_4gBw

% 200 ns
% 2 GHz = .5 ns  - 500ps, 400 pulses in 200ns
% 5e3pts, 40ps/pt, 12.5 pts per pulse,  
% interp - 10e3 pts, 20ps/pt, 25pts per pulse

% 3 GHz = .333ns, 600 pulses in 200 ns
% 5e3pts, 40ps/pt, 8.33 pts per pulse,

% 4 GHz = .25ns, 800 pulses in 200 ns
% 5e3pts, 40ps/pt, 6.25 pts per pulse,


% Calculating the One Level
% Middle 20% of data
% 50ps to 75ps or -25ps to 25ps

% Get middle 20% of data
low = find(time_eye>-31e-12, 1);
high = find(time_eye<31e-12, 1, 'last');
middledata = eye(low:high, :);

% Seperate data into one and zero levels
one_lvl = middledata(:,middledata(1,:) > 0);
one_lvl_avg = mean(one_lvl);

zero_lvl = middledata(:,middledata(1,:) < 0);
zero_lvl_avg = mean(zero_lvl);

% Calculate I1 mean and std
one_lvl_mean = mean(one_lvl_avg);
one_lvl_std = std(one_lvl_avg);


h = figure();
set(h,'WindowStyle','docked');

plothist(one_lvl_avg, 10, mean(one_lvl_avg), std(one_lvl_avg))
title("One-Level Voltage Distribution")
xlabel("One-Level Voltage [V]")
ylabel("PDF")
legend("Experimental Data", "Gaussian Curve Fit")
% print -depsc onelvl_2gSig_4gBw


% histogram(one_lvl_avg, 10, 'Normalization','pdf'); hold on;
% 
% y = min(one_lvl_avg):0.0001:max(one_lvl_avg);
% mu = mean(one_lvl_avg);
% sigma = std(one_lvl_avg);
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5);
% xline(mean(one_lvl_avg))

% Calculate I0 mean and std
zero_lvl_mean = mean(zero_lvl_avg);
zero_lvl_std = std(zero_lvl_avg);

h = figure();
set(h,'WindowStyle','docked');

% histogram(zero_lvl_avg, 10, 'Normalization','pdf'); hold on;
% 
% y = -.265:0.0001:-.22;
% mu = zero_lvl_mean;
% sigma = zero_lvl_std;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)
% title("Zero Level Distribution")
% xlabel("Zero Level Voltage Level [V]")
% ylabel("PDF")

plothist(zero_lvl_avg, 10, zero_lvl_mean, zero_lvl_std)
title("Zero-Level Distribution")
xlabel("Zero-Level Voltage Level [V]")
ylabel("PDF")
legend("Experimental Data", "Gaussian Curve Fit")
% print -depsc zerolvl_2gSig_4gBw

% Calculate Eye Amplitude
eye_amp = one_lvl_mean - zero_lvl_mean;
eye_ht = one_lvl_mean - 3*one_lvl_std - zero_lvl_mean - 3*zero_lvl_std;

% Calculate Q value
Q = eye_amp/(one_lvl_std+zero_lvl_std);

% First Crossing points
rise = 0;
crossT = findcross(time_eye, eye, 1, 0, rise);

% Calculate CrossT mean and std
crossT_mean = mean(crossT);
crossT_std = std(crossT);

jitter_rms = crossT_std;
jitter_pp = max(crossT)-min(crossT);

h = figure();
set(h,'WindowStyle','docked');

% histogram(crossT, 10, 'Normalization','pdf'); hold on;
% 
% y = -1.26e-10:1e-13:-1.1e-10;
% mu = crossT_mean;
% sigma = crossT_std;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5);

plothist(crossT, 10, crossT_mean, crossT_std)
title("Jitter Distribution");
xlabel("Time of 50% Crossing [s]");
ylabel("PDF");
legend("Experimental Data", "Gaussian Curve Fit")
% print -depsc jitter_2gSig_4gBw

% Second Crossing points
crossT2 = findcross(time_eye, eye, 0, 0, rise);

% Calculate CrossT2 mean and std
crossT2_mean = mean(crossT2);
crossT2_std = std(crossT2);

h = figure();
set(h,'WindowStyle','docked');

% histogram(crossT2, 10, 'Normalization','pdf'); hold on;
% 
% y = 1.26e-10:1e-13:1.36e-10;
% mu = crossT2_mean;
% sigma = crossT2_std;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5);

plothist(crossT2, 10, crossT2_mean, crossT2_std)
title("Jitter Distribution");
xlabel("Time of 50% Crossing [s]");
ylabel("PDF");

legend("Experimental Data", "Gaussian Curve Fit")

% Eye Width
eye_width = (crossT2_mean - 3*crossT2_std) - (crossT_mean + 3*crossT_std);

% % Rise time
% % Calculate 10 and 90% Level
% pct10 = 0.1*eye_amp + zero_lvl_mean;
% pct90 = 0.9*eye_amp + zero_lvl_mean;
% 
% rise = 1;
% rise_pct10_cross = findcross(time_eye, eye, 1, pct10, rise);
% rise_pct10_cross_mean = mean(rise_pct10_cross);
% rise_pct10_cross_std = std(rise_pct10_cross);
% 
% rise_pct90_cross = findcross(time_eye, eye, 1, pct90, rise);
% rise_pct90_cross_mean = mean(rise_pct90_cross);
% rise_pct90_cross_std = std(rise_pct90_cross);
% 
% h = figure();
% set(h,'WindowStyle','docked');
% 
% histogram(rise_pct10_cross, 10, 'Normalization','pdf'); hold on;
% histogram(rise_pct90_cross, 10, 'Normalization','pdf'); hold on;
% 
% % y = -1.7e-10:1e-13:-1.4e-10;
% mu = rise_pct10_cross_mean;
% sigma = rise_pct10_cross_std;
% y = (mu-4*sigma):abs(mu)/1e4:(mu+4*sigma);
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',2,'color', 'b');
% 
% % y = -.9e-10:1e-13:-.7e-10;
% mu = rise_pct90_cross_mean;
% sigma = rise_pct90_cross_std;
% y = (mu-4*sigma):abs(mu)/1e4:(mu+4*sigma);
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',2,'color', 'r');
% 
% title("Rise Time Distribution");
% xlabel("Time of 10% and 90% Crossing [s]");
% ylabel("PDF");
% legend("10% Crossing", "90% Crossing", 'location', 'north')
% 
% % print -depsc risetime_2gSig_4gBw
% 
% risetime = rise_pct90_cross_mean - rise_pct10_cross_mean;
% 
% rise = 2;
% fall_pct10_cross = findcross(time_eye, eye, 1, pct90, rise);
% fall_pct10_cross_mean = mean(fall_pct10_cross);
% fall_pct10_cross_std = std(fall_pct10_cross);
% 
% fall_pct90_cross = findcross(time_eye, eye, 1, pct10, rise);
% fall_pct90_cross_mean = mean(fall_pct90_cross);
% fall_pct90_cross_std = std(fall_pct90_cross);
% 
% h = figure();
% set(h,'WindowStyle','docked');
% 
% histogram(fall_pct10_cross, 10, 'Normalization','pdf'); hold on;
% histogram(fall_pct90_cross, 10, 'Normalization','pdf'); hold on;
% 
% mu = fall_pct10_cross_mean;
% sigma = fall_pct10_cross_std;
% y = (mu-4*sigma):abs(mu)/1e4:(mu+4*sigma);
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',2,'color', 'b');
% 
% mu = fall_pct90_cross_mean;
% sigma = fall_pct90_cross_std;
% y = (mu-4*sigma):abs(mu)/1e4:(mu+4*sigma);
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',2,'color', 'r');
% 
% xlim([(fall_pct10_cross_mean-4*fall_pct10_cross_std), (mu+4*sigma)])
% 
% title("Fall Time Distribution");
% xlabel("Time of 10% and 90% Crossing [s]");
% ylabel("PDF");
% legend("10% Crossing", "90% Crossing", 'location', 'north')
% 
% % print -depsc falltime_2gSig_4gBw
% 
% falltime = fall_pct90_cross_mean - fall_pct10_cross_mean;

names = {'Zero Level mu', 'Zero Level std', 'One Level mu', 'One Level std', ...
    'Eye amp', 'Eye ht', 'eye width', 'Q', 'jitter rms', 'jitter pp'};
%'rise time', 'fall time',
% risetime, falltime,
values = [zero_lvl_mean, zero_lvl_std, one_lvl_mean, one_lvl_std, eye_amp, ...
    eye_ht, eye_width, Q, jitter_rms, jitter_pp];

format long

values
% for i=1:length(names)
%      fprintf('%s=%f\n',names{i},values(i))
% end