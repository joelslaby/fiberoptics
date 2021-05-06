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

filename = '\media\lab5_day3_000_ALL.csv';
    
T = readtable(filename);
data = table2array(T);

time = data(6:end, 1);
time_check = time;
wave = data(6:end, 2);
interp_factor = 1;
time_interp = linspace(time(1), time(end), interp_factor*length(time)-1);
wave = interp1(time, wave, time_interp, 'spline');
time = time_interp;
% wave_filt = lowpass(wave, 1e9, 25e9);

%Clock Recovery
%Grabs max and min and creates thresholds to find transitions
wave_max = max(wave);
wave_min = min(wave);
th_1 = (wave_max+wave_min)/2;
th_0 = (th_1 + wave_min)/2;
th_2 = (th_1 + wave_max)/2;
%Try different threshold creation methods
% th_1 = (wave_max+wave_min)/2;
% th_0 = wave_min/2;
% th_2 = wave_max/2;
%Loop 'discretizes' values and ensures transitions aren't falsly made into
%a descretized point. 
prev_state = 0;
cur_state = 0;
state_ctr = 10;
for i = 1:length(wave)
%     if wave(i) > th_2
%         
%         wave_th(i) = 0.1;
%         state = 1;
%     elseif (wave(i)<th_2 && wave(i)>th_1)
%         if(state ~=2)
%             if(state_ctr > 2)
%                 Normal transition
%                 wave_th(i) = 0.05;
%             else
%                 Found false transition
%                 
%             end
%             state_ctr = 0;
%         else
%             Same state
%             wave_th(i) = 0.05;
%             state_ctr = state_ctr + 1;
%         end
%         state = 2;
%     elseif (wave(i)<th_1 && wave(i)>th_0)
%         cur_state = th_0;
%         if(prev_state ~= cur_state)
%             wave_th(i) = cur_state;
%             
%             if(state_ctr <= 2)
%                 Found false transition
%                 wave_th(i-2) = wave_th(i-3);
%                 wave_th(i-1) = cur_state;
%             end
%             state_ctr = 0;
%         else
%             Same state
%             wave_th(i) = 0.05;
%             state_ctr = state_ctr + 1;
%         end
%         prev_state = 2;
%         
%         
%         wave_th(i) = -0.05;
%     else
%         wave_th(i) = -0.1;
%     end
    
    if wave(i) > th_2
        cur_state = 3;
    elseif (wave(i)<th_2 && wave(i)>th_1)
        cur_state = 2;
    elseif (wave(i)<th_1 && wave(i)>th_0)
        cur_state = 1;
    else
        cur_state = 0;
    end
    
    if(i < 4*interp_factor)
        prev_state = cur_state;
    end
    
    if(prev_state ~= cur_state)
        wave_th(i) = cur_state;

        if(state_ctr < 2*interp_factor)
            % Found false transition
            wave_th(i-3*interp_factor:i-1*interp_factor-1) = wave_th(i-(3*interp_factor+1));
            wave_th(i-1*interp_factor:i-1) = cur_state;
%             wave_th(i-3:i-2) = wave_th(i-4);
%             wave_th(i-1) = cur_state;
        end
        state_ctr = 0;
    else
        % Same state
        wave_th(i) = cur_state;
        state_ctr = state_ctr + 1;
    end
    prev_state = cur_state; 
end

%Converts descretized values into marking the transitions and finally
%finding the time difference between each transition, thus giving us a list
%of bit periods to process.
time_changes = conv(wave_th, [1 -1]);
time_changes(time_changes ~= 0) = 1;
time_changes = logical(time_changes(1:end-1));
time_differences = conv(time(time_changes), [1 -1]);
time_differences = time_differences(2:end-1);

%Oldway (Wrapped in comments)
wave_th1 = zeros(1,length(wave));
wave_th1(wave>th_2) = 0.06;
wave_th1(wave<th_2 & wave>th_1) = 0.04;
wave_th1(wave<th_1 & wave>th_0) = 0.02;
wave_th1(wave<th_0) = 0;
time_changes1 = conv(wave_th1, [1 -1]);
time_changes1(1 ~= 0) = 1;
time_changes1 = logical(time_changes1(1:end-1));
time_differences1 = conv(time(time_changes1), [1 -1]);
time_differences1 = time_differences1(2:end-1);
% time_differences = time_differences/max(time_differences);
%End Old way (End wrapper)

time_cutoff = 0.00; %Trims min bit periods found to remove errors
% time_differences = time_differences(time_differences>(max(time_differences)*time_cutoff));
time_differences = time_differences(time_differences>time_cutoff*3*41e-12); %Can cut out super low unexpected results
time_differences = sort(time_differences);

time_differences_norm = time_differences/min(time_differences(floor(0.5*length(time_differences)):end));

bin = time_differences(time_differences_norm > 0.5 & time_differences_norm < 1.5);
freq_recov = 1/mean(bin)
histogram(time_differences/1e-9, 500);
title("Distribution of Time Differences");
xlabel("Time Differences [ns]");
ylabel("Number of Occurences")

%Plot
ps = 1e-12;
mV = 1e-3;
h = figure();
set(h,'WindowStyle','docked');
% subplot(2, 1, 1);
Ts = 40e-12;
endNbr = 200;
ns = 1e-9;
time_sample = (1:endNbr)*Ts;

subplot(2, 1, 1);
% plot(time(1:endNbr)/ps, wave(1:endNbr)/mV,time(1:endNbr)/ps, wave_filt(1:endNbr)/mV); hold on;
plot(time_sample/ns, wave(15:14+endNbr)/mV, 'linewidth', 2);hold on;
xlim([-1, 8]);
yline(th_0/mV);
yline(th_1/mV);
yline(th_2/mV);
subplot(2, 1, 2);
plot(time_sample/ns, (wave_th(15:14+endNbr)), 'linewidth', 2);
% plot(time_sample/ns, wave_th1(1:endNbr)/mV); 
% plot(time(1:endNbr)/ps, wave(1:endNbr)/mV,time(1:endNbr)/ps, wave_th(1:endNbr)/mV); hold on;
% plot(time(1:endNbr)/ps, wave(1:endNbr)/mV,time(1:endNbr)/ps, wave_th(1:endNbr)/mV,time(1:endNbr)/ps, time_changes(1:endNbr)/mV); hold on;
xlim([-1, 8]);
yline(0.5);
yline(1.5);
yline(2.5);
% legend('original','corrected','shitty', 'Thresh 1', 'Thresh 2', 'Thresh 3');
% yline(th_3/mV);
diff = [-1, 1];
ylim([-.5, 3.5])

waveDiff = conv(wave, diff);

% subplot(2, 1, 2);
% plot(time(1:endNbr)/ps, waveDiff(1:endNbr)/mV);

% pk = max(waveDiff);
% waveDiff(waveDiff>.7*pk) = 1;
% waveBin = waveDiff > .6*max(waveDiff) =
% plot(time(1:endNbr)/ps, waveDiff(1:endNbr)/mV);
% findpeaks(waveDiff(1:endNbr)/mV);
% h = figure();
% set(h,'WindowStyle','docked');
% % plot(time, )
% hold on
% yline(0.1)
h = figure();
set(h,'WindowStyle','docked');
% plot(time_differences); hold on;
plot(time_differences_norm);

%%

h = figure();
set(h,'WindowStyle','docked');
Ts = 40e-12;
endNbr = 200;
ns = 1e-9;
time_sample = (1:endNbr)*Ts;

plot(time_sample/ns, wave(15:14+endNbr)/mV,'k', 'linewidth', 2); hold on;
% xlim([time(1)/ps, time(endNbr)/ps]);
yline(th_0/mV, 'r','linewidth', 2);
yline(th_1/mV, 'g','linewidth', 2);
yline(th_2/mV, 'b','linewidth', 2);
plot(time_sample/ns, 100*time_changes(15:14+endNbr)-250, 'linewidth', 2)
% legend('original','corrected','Thresh 1', 'Thresh 2', 'Thresh 3');
ylim([-250, 250]);
xlim([-1, 8]);
text(-.5, -150, "0",'FontSize',14);
text(-.5, -50, "1",'FontSize',14);
text(-.5, 50, "2",'FontSize',14);
text(-.5, 150, "3",'FontSize',14);
xlabel("Time [ns]");
ylabel("Amplitude [mV]");
title("Received Signal with Thresholds");

% print -depsc received_sig_thresholds
% exportgraphics(h,'received_sig_thresholds.png','Resolution',300);

