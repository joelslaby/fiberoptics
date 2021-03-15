clearvars; clc;
close all;

% Import waveform data
T = readtable('media/lab4_007_ALL.csv');
data = table2array(T);

% Separate data to time and wave vectors
time = data(6:end, 1);
wave = data(6:end, 3);

% Signal Information
setup = 0;
showplots = ~setup;
% 2 GHz = .5 ns  - 500ps, 4000 pulses in 2us
% 50e3pts, 40ps/pt, 12.5 pts per pulse,
% interp - 400e3 pts, 20ps/pt, 25pts per pulse
sig_freq = 4; %GHz 
clk_freq = 2*sig_freq;
% Lab4 - 000->40
%        001->45
%        002->43
%        003->10
%        004->8
%        005->18
%        006->30
%        007->25
%        008->43
offset = 25;
packet = 50;
nbr_traces = 1998;
nbr_bins = 100;

% Scaling factors
ps = 1e-12;
mV = 1e-3;

% Interpolate Data
time_interp = linspace(time(1), time(end), length(time)*(2*clk_freq)-1);
wave_interp = interp1(time, wave, time_interp, 'spline');

% Create time vector
time_eye = ((1:packet*2)-packet)*20e-12/clk_freq;
time_eye_ps = time_eye/ps;

% Create eye matrix
eye = zeros(packet*2, nbr_traces);
for i=1:nbr_traces
    eye(:, i) = wave_interp((i*packet+offset-packet/2):(((i+1)*packet+offset+packet/2)-1))';
end

eye_mv = eye/mV;

% Create Eye Gradient matrix    
eye_grad = zeros(2*packet, nbr_bins);
edges = linspace(min(min(eye)), max(max(eye)), nbr_bins+1);
bins = edges(1:end-1) + (edges(2)-edges(1))/2;
[X, Y] = meshgrid(time_eye_ps, bins);
for i=1:2*packet
    [N,edge,bin] = histcounts(eye(i, :), edges);
    eye_grad(i, :) = N';
end

% Log scale eye gradient
scaled_eye = 10*log10(eye_grad);
scaled_eye(scaled_eye == -Inf) = 0;

if(setup)
%     h = figure();
%     set(h,'WindowStyle','docked');
    surf(X, Y,scaled_eye');
    view(0,90);
    shading interp;
    xlim([time_eye_ps(1), time_eye_ps(end)]);
    cmap = colormap(turbo);
    colorbar;
    set(gca,'Color',cmap(1, :))
end

% Calculating the One and Zero Level
% Middle 20% of data
ten_pct = 1 / (sig_freq * 1e9) / 2 * .12;
low = find(time_eye>-ten_pct, 1);
high = find(time_eye<ten_pct, 1, 'last');
middledata = eye(low:high, :);

% Seperate data into one and zero levels
midpt = .015;
one_lvl = middledata(:,middledata(1,:) > midpt);
one_lvl_avg = mean(one_lvl);

zero_lvl = middledata(:,middledata(1,:) < midpt);
zero_lvl_avg = mean(zero_lvl);

% Calculate I1 mean and std
one_lvl_mean = mean(one_lvl_avg);
one_lvl_std = std(one_lvl_avg);

% Calculate I0 mean and std
zero_lvl_mean = mean(zero_lvl_avg);
zero_lvl_std = std(zero_lvl_avg);

% Calculate Eye Amplitude
eye_amp = one_lvl_mean - zero_lvl_mean;
eye_ht = one_lvl_mean - 3*one_lvl_std - zero_lvl_mean - 3*zero_lvl_std;

% Calculate Q value
Q = eye_amp/(one_lvl_std+zero_lvl_std);

% Find waveform crossing points
pct20 = 0.2*eye_amp + zero_lvl_mean;
pct50 = (one_lvl_mean+zero_lvl_mean)/2
pct80 = 0.8*eye_amp + zero_lvl_mean;
[crossing middlepts] = findcross(time_eye, eye, packet, [pct20, pct50, pct80]);
L = length(crossing);

pct50_1 = InterX([time_eye_ps, eye(:, middlepts(1, 1))';
    time_eye_ps, eye(:, middlepts(3, 1))']);
pct50_2 = InterX([time_eye_ps, eye(:, middlepts(2, 1))';
    time_eye_ps, eye(:, middlepts(4, 1))']);
pct50 = mean([pct50_1(1), pct50_2(1)]);

crossing_pct = (pct50 - zero_lvl_mean)/(one_lvl_mean - zero_lvl_mean);

[crossing middlepts] = findcross(time_eye, eye, packet, [pct20, pct50, pct80]);
L = length(crossing);

% First Crossing points
size(crossing(1, 1:2, 2, :))
L
crossT = reshape(crossing(1, 1:2, 2, :), [1, 2*L]);

% Calculate CrossT mean and std
crossT_mean = mean(crossT);
crossT_std = std(crossT);

jitter_rms = crossT_std;
jitter_pp = max(crossT)-min(crossT);

% Second Crossing points
crossT2 = reshape(crossing(2, 1:2, 2, :), [1, 2*L]);

% Calculate CrossT2 mean and std
crossT2_mean = mean(crossT2);
crossT2_std = std(crossT2);

% Eye Width
eye_width = (crossT2_mean - 3*crossT2_std) - (crossT_mean + 3*crossT_std);

% Rise time
rise_pct20_cross = reshape(crossing(1, 1, 1, :), [1, L]);
rise_pct80_cross = reshape(crossing(1, 1, 3, :), [1, L]);

rise_time = rise_pct80_cross - rise_pct20_cross;
rise_time_mean = mean(rise_time);
rise_time_std = std(rise_time);

% Fall time
fall_pct20_cross = reshape(crossing(1, 2, 1, :), [1, L]);
fall_pct80_cross = reshape(crossing(1, 2, 3, :), [1, L]);
fall_time = fall_pct20_cross - fall_pct80_cross;
fall_time_mean = mean(fall_time);
fall_time_std = std(fall_time);

% Plots
if(showplots)
    % Plot time waveform
    h = figure();
    set(h,'WindowStyle','docked');
    plot(time, wave); hold on;
    plot(time_interp, wave_interp);
    title("Eye Diagram - 6GHz signal - 4 GHz Bandwidth");
    xlabel("Time [ps]");
    ylabel("Voltage [V]");

    % Plot eye diagram
    h = figure();
    set(h,'WindowStyle','docked');
    plot(time_eye_ps, eye_mv);
    xlim([time_eye_ps(1), time_eye_ps(end)]);
    title("Eye Diagram - 6GHz signal - 4 GHz Bandwidth");
    xlabel("Time [ps]");
    ylabel("Voltage [V]");

    % Plot eye diagram gradient
    h = figure();
    set(h,'WindowStyle','docked');
    surf(X, Y,scaled_eye');
    view(0,90);
    shading interp;
    xlim([time_eye_ps(1), time_eye_ps(end)]);
    cmap = colormap(turbo);
    colorbar;
    set(gca,'Color',cmap(1, :))
    
    % Plot one-level voltage distribution
    h = figure();
    set(h,'WindowStyle','docked');
    plothist(one_lvl_avg, 10)
    title("One-Level Voltage Distribution")
    xlabel("One-Level Voltage [V]")
    ylabel("PDF")
    legend("Experimental Data", "Gaussian Curve Fit")

    % Plot zero-level voltage distribution
    h = figure();
    set(h,'WindowStyle','docked');
    plothist(zero_lvl_avg, 10)
    title("Zero-Level Voltage Distribution")
    xlabel("Zero-Level Voltage Level [V]")
    ylabel("PDF")
    legend("Experimental Data", "Gaussian Curve Fit")

    % Plot Jitter distribution
    h = figure();
    set(h,'WindowStyle','docked');
    plothist(crossT, 10)
    title("Jitter Distribution");
    xlabel("Time of 50% Crossing [s]");
    ylabel("PDF");
    legend("Experimental Data", "Gaussian Curve Fit")
    
    % Plot Rise Time distribution
    h = figure();
    set(h,'WindowStyle','docked');
    plothist(rise_time, 10)
    title("Rise Time Distribution");
    xlabel("Rise Time [ps]");
    ylabel("PDF");
    
    % Plot Fall Time distribution
    h = figure();
    set(h,'WindowStyle','docked');
    plothist(fall_time, 10)
    title("Fall Time Distribution");
    xlabel("Fall Time [ps]");
    ylabel("PDF");
end
