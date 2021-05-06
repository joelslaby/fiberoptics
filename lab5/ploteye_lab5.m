clearvars; clc;
close all;

%% CODE
clearvars; clc;
close all;
% SETUP
% 0 = No plots, Calculate param 
% 1 = Plot Gradient Eye Only, No param calculations, Used to center graph
% 2 = Plot everything, Calculate all param
% 3 = Plot all traces in Grid
plot_data = 2;

% -1-select fiber
% 0-no fiber, no mod
% 1-no fiber, mod
% 2-fiber, mod
% 3-lab4 day1, wrong vPI
% 4 - lab5, day4
wave_type = 4;
select = 0;

% Initialize Data Arrays
Q_save = zeros(1, 21);
BER_save = zeros(1, 21);
one_lvl_mean_save = zeros(1, 21);
zero_lvl_mean_save = zeros(1, 21);
mid0_lvl_mean_save = zeros(1, 21);
mid1_lvl_mean_save = zeros(1, 21);
one_lvl_std_save = zeros(1, 21);
zero_lvl_std_save = zeros(1, 21); 
mid0_lvl_std_save = zeros(1, 21);
mid1_lvl_std_save = zeros(1, 21); 

one_lvl_mean_norm_save = zeros(1, 21);
zero_lvl_mean_norm_save = zeros(1, 21);
mid0_lvl_mean_norm_save = zeros(1, 21);
mid1_lvl_mean_norm_save = zeros(1, 21);
one_lvl_std_norm_save = zeros(1, 21);
zero_lvl_std_norm_save = zeros(1, 21); 
mid0_lvl_std_norm_save = zeros(1, 21);
mid1_lvl_std_norm_save = zeros(1, 21);

freq_recov_save = zeros(1, 21);

Q_01_save = zeros(1, 21);
Q_02_save = zeros(1, 21);
Q_03_save = zeros(1, 21);
Q_12_save = zeros(1, 21);
Q_13_save = zeros(1, 21);
Q_23_save = zeros(1, 21);

if(select)
    nbrTraces = 1:1;
elseif(wave_type > 0)
    nbrTraces = 1:21;
else
    nbrTraces = 1:8;
end

if(plot_data == 3)
    h = figure();
    h.Position = [100 100 1500 900];
%     set(h,'WindowStyle','docked');
end

a = 1;
b = 1;
c = 1;

for p=nbrTraces
    % Import waveform data
    if(wave_type ~= -1)
        if(wave_type < 2)
            selectdata = num2str(p+8, '%.2d');
        elseif(wave_type == 2)
            selectdata = num2str(p+17, '%.2d');
        else
            selectdata = num2str(p-1, '%.2d');
        end
        if(wave_type < 4) 
            filename = strcat('media_day2/lab4_0',selectdata,'_ALL.csv');
        else
            filename = strcat('media\lab5_day4_0',selectdata,'_ALL.csv');
        end
    else
        filename = '\media\lab5_day4_020_ALL.csv';
    end
    T = readtable(filename);
    data = table2array(T);

    % Separate data to time and wave vectors
    time = data(6:end, 1);
    if(wave_type < 0 & wave_type < 4)
        wave = data(6:end, 3);
    else 
        wave = data(6:end, 2);
    end

    % 2 GHz = .5 ns  - 500ps, 4000 pulses in 2us
    % 50e3pts, 40ps/pt, 12.5 pts per pulse,
    % interp - 400e3 pts, 20ps/pt, 25pts per pulse
    if(wave_type == -1)
        sig_freq = 3.9;
    elseif(wave_type < 2)
        sig_freqs = [2, 2, 2, 3, 3, 3, 4, 4, 4];
        sig_freq = sig_freqs(p); %GHz
    elseif(wave_type == 4)
        sig_freqs = [0.975, 1.95, 3.9];
        sig_freqs = repmat(sig_freqs, [1, 7]);
        sig_freq = sig_freqs(p); %GHz        
    else
        sig_freq = 2;
    end
    
    clk_freq = sig_freq;
    
    % Determine offsets to center graph
    if(wave_type == 0)
        sig_offsets = [25, 27, 27, 15, 5, 3, 36, 2]; % Lab 4, no fiber, no mod
    elseif(wave_type == 1)
        sig_offsets = [6, 10, 8, 40, 29, 28, 0, 17, 47]; % Lab 4, no fiber
    elseif(wave_type == 2)
        sig_offsets = [46, 38, 18, 32, 10, 20, 31, 46, 35]; % Lab 4 w/ fiber
    elseif(wave_type == 4)
        sig_offsets = [0,5,28,35,8,0,45,42,23,33,2,12,46,40,30,46,40,46,46,41,29];
    else
        sig_offsets = 0;
    end
    
    offset = sig_offsets(p);
    
    % Waveform properties
    packet = 50;
    nbr_traces = 999*round(sig_freq, 0);
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
    edges = linspace(min(min(eye_mv)), max(max(eye_mv)), nbr_bins+1);
    bins = edges(1:end-1) + (edges(2)-edges(1))/2;
    [X, Y] = meshgrid(time_eye_ps, bins);
    for i=1:2*packet
        [N,edge,bin] = histcounts(eye_mv(i, :), edges);
        eye_grad(i, :) = N';
    end
    
    % Log scale eye gradient
    scaled_eye = 10*log10(eye_grad);
    scaled_eye(scaled_eye == -Inf) = 0;

    if((plot_data ~= 1) & (plot_data ~= 3))
        % Calculating the One and Zero Level
        % Middle 24% of data
        ten_pct = 1 / (sig_freq * 1e9) / 2 * .12;
        low = find(time_eye>-ten_pct, 1);
        high = find(time_eye<ten_pct, 1, 'last');
        middledata = eye(low:high, :);

        % Seperate data into one and zero levels
        min_eye = min(min(eye));
        max_eye = max(max(eye));
        midpt = (min_eye + max_eye) / 2;
%         high_thresh = (midpt + 1.1 * max_eye)/2
        high_thresh = .4*midpt + .6*max_eye;
        one_lvl = middledata(:,middledata(1,:) > high_thresh);
        one_lvl_avg = mean(one_lvl);

%         low_thresh = (midpt + 0.5 * min_eye)/(1+.5)
        low_thresh = .4*midpt + .6*min_eye;
        zero_lvl = middledata(:,middledata(1,:) < low_thresh);
        zero_lvl_avg = mean(zero_lvl);
        
        mid0_lvl = middledata(:,(middledata(1,:) > low_thresh) & (middledata(1,:) < midpt));
        mid0_lvl_avg = mean(mid0_lvl);
        
        mid1_lvl = middledata(:,(middledata(1,:) < high_thresh) & (middledata(1,:) > midpt));
        mid1_lvl_avg = mean(mid1_lvl);

        % Calculate I1 mean and std
        one_lvl_mean = mean(one_lvl_avg);
        one_lvl_std = std(one_lvl_avg);

        one_lvl_mean_norm = mean(one_lvl_avg/one_lvl_mean);
        one_lvl_std_norm = std(one_lvl_avg/one_lvl_mean);
        
        % Calculate I0 mean and std
        zero_lvl_mean = mean(zero_lvl_avg);
        zero_lvl_std = std(zero_lvl_avg);   
        
        zero_lvl_mean_norm = mean(zero_lvl_avg/one_lvl_mean);
        zero_lvl_std_norm = std(zero_lvl_avg/one_lvl_mean);
 
        mid0_lvl_mean = mean(mid0_lvl_avg);
        mid0_lvl_std = std(mid0_lvl_avg);

        mid0_lvl_mean_norm = mean(mid0_lvl_avg/one_lvl_mean);
        mid0_lvl_std_norm = std(mid0_lvl_avg/one_lvl_mean);
        
        mid1_lvl_mean = mean(mid1_lvl_avg);
        mid1_lvl_std = std(mid1_lvl_avg);
        
        mid1_lvl_mean_norm = mean(mid1_lvl_avg/one_lvl_mean);
        mid1_lvl_std_norm = std(mid1_lvl_avg/one_lvl_mean);
        
        % Save Data
        one_lvl_mean_save(p) = one_lvl_mean;
        zero_lvl_mean_save(p) = zero_lvl_mean;
        mid0_lvl_mean_save(p) = mid0_lvl_mean;
        mid1_lvl_mean_save(p) = mid1_lvl_mean;
        one_lvl_std_save(p) = one_lvl_std;
        zero_lvl_std_save(p) = zero_lvl_std; 
        mid0_lvl_std_save(p) = mid0_lvl_std;
        mid1_lvl_std_save(p) = mid1_lvl_std; 
        
        one_lvl_mean_norm_save(p) = one_lvl_mean_norm;
        zero_lvl_mean_norm_save(p) = zero_lvl_mean;
        mid0_lvl_mean_norm_save(p) = mid0_lvl_mean_norm;
        mid1_lvl_mean_norm_save(p) = mid1_lvl_mean_norm;
        one_lvl_std_norm_save(p) = one_lvl_std_norm;
        zero_lvl_std_norm_save(p) = zero_lvl_std_norm; 
        mid0_lvl_std_norm_save(p) = mid0_lvl_std_norm;
        mid1_lvl_std_norm_save(p) = mid1_lvl_std_norm; 

        % Calculate Eye Amplitude
        eye_amp = one_lvl_mean - zero_lvl_mean;
        eye_ht = one_lvl_mean - 3*one_lvl_std - zero_lvl_mean - 3*zero_lvl_std;

        eye_amp_01 = mid0_lvl_mean - zero_lvl_mean;
        eye_amp_02 = mid1_lvl_mean - zero_lvl_mean;
        eye_amp_03 = one_lvl_mean - zero_lvl_mean;
        eye_amp_12 = mid1_lvl_mean - mid0_lvl_mean;
        eye_amp_13 = one_lvl_mean - mid0_lvl_mean;
        eye_amp_23 = one_lvl_mean - mid1_lvl_mean;
        
        Q_01 = eye_amp_01/(mid0_lvl_std + zero_lvl_std);
        Q_02 = eye_amp_02/(mid1_lvl_std + zero_lvl_std);
        Q_03 = eye_amp_03/(one_lvl_std + zero_lvl_std);
        Q_12 = eye_amp_12/(mid1_lvl_std + mid0_lvl_std);
        Q_13 = eye_amp_13/(one_lvl_std + mid0_lvl_std);
        Q_23 = eye_amp_23/(one_lvl_std + mid1_lvl_std);
        
        Q_01_save(p) = Q_01';
        Q_02_save(p) = Q_02';
        Q_03_save(p) = Q_03';
        Q_12_save(p) = Q_12';
        Q_13_save(p) = Q_13';
        Q_23_save(p) = Q_23';
        
        % Calculate Q value and BER
%         Q = eye_amp/(one_lvl_std+zero_lvl_std);
%         BER = 1/2 * erfc(Q./sqrt(2));
        
%         Q_save(p) = Q';
%         BER_save(p) = BER';

        %Find New wave crossing points
        midpt = (mid0_lvl_std*mid1_lvl_mean + mid1_lvl_std*mid0_lvl_mean)...
            /(mid0_lvl_std+mid1_lvl_std);
        high_thresh = (mid1_lvl_std*one_lvl_mean + one_lvl_std*mid1_lvl_mean)...
            /(one_lvl_std+mid1_lvl_std);
        low_thresh = (mid0_lvl_std*zero_lvl_mean + zero_lvl_std*mid0_lvl_mean)...
            /(zero_lvl_std+mid0_lvl_std);
%  
%         midpt = (mid1_lvl_mean + mid0_lvl_mean)...
%             /(2);
%         high_thresh = (one_lvl_mean + mid1_lvl_mean)...
%             /(2);
%         low_thresh = (zero_lvl_mean + mid0_lvl_mean)...
%             /(2);
        % Recover Bit Stream
        avg_middle = mean(middledata);
        bit_stream = zeros(1, nbr_traces);
        bit_stream(find(avg_middle > high_thresh)) = 1;
        bit_stream(find((avg_middle < high_thresh) & (avg_middle > midpt))) = 1/3;
        bit_stream(find((avg_middle < midpt)&(avg_middle > low_thresh))) = -1/3;
        bit_stream(find(avg_middle < low_thresh)) = -1;
        
        nbrBits = length(bit_stream);
        fDrive = sig_freq;
        tPulse = 1/sig_freq;
        ptsPerPulse = 50;

        % Adds pts at 0.5 so that it transitions perfectly on each bit
        diff = abs(conv(bit_stream, [1 -1]));
        trig = reshape([diff; zeros(ptsPerPulse-1, length(diff))], 1, []);
        
        % Create the drive signal time vector
        tDrive = linspace(0, nbrBits * tPulse, nbrBits * ptsPerPulse + 1);

        % Create the drive signal
        driveSig = zeros(1, ptsPerPulse*nbrBits + 1);
        driveSig(1:end-1) = repelem(bit_stream, ptsPerPulse);
        driveSig(find(trig == 1)) = .5;
        driveSig = (one_lvl_mean-zero_lvl_mean)/2*driveSig+(one_lvl_mean+zero_lvl_mean)/2;
        
        % Create the PRBS Signal
        O = 15;
        offset = 3*round(sig_freq, 0);

        L = 2^(O+1)-1;

        sig = prbs(O,L);
        sig = sig*2-1;
        sigP = -0.5 * circshift(sig, -offset);
        sumSig = (sig+sigP);

        bit_stream_prbs = zeros(1, length(sig));
        bit_stream_prbs(find(sumSig == -1.5)) = -1;
        bit_stream_prbs(find(sumSig == -.5)) = 1/3;
        bit_stream_prbs(find(sumSig == .5)) = -1/3;
        bit_stream_prbs(find(sumSig == 1.5)) = 1;
        
        bit_corr = conv(bit_stream_prbs, flip(bit_stream));
        [corr_max prbs_start] = max(bit_corr);
       
        
        shifted_prbs = circshift(bit_stream_prbs, -(prbs_start - length(bit_stream)));
        shifted_prbs = shifted_prbs(1:length(bit_stream));
%         shifted_prbs = bit_stream_prbs((length(bit_corr)-prbs_start):(length(bit_corr)-prbs_start+length(bit_stream)-1));
        
        prbsSig = zeros(1, ptsPerPulse*nbrBits + 1);
        prbsSig(1:end-1) = repelem(shifted_prbs, ptsPerPulse);
        prbsSig(find(trig == 1)) = .5;
        prbsSig = (one_lvl_mean-zero_lvl_mean)/2*prbsSig+(one_lvl_mean+zero_lvl_mean)/2;
        
        find((shifted_prbs - bit_stream) ~= 0);
        BER = length(find((shifted_prbs - bit_stream) ~= 0))/length(bit_stream);
        BER_save(p) = BER;
        
        %Clock Recovery
        %Grabs max and min and creates thresholds to find transitions
%         wave_max = max(wave);
%         wave_min = min(wave);
%         th_1 = (wave_max+wave_min)/2; %midpt
%         th_0 = (th_1 + wave_min)/2; % low_thresh
%         th_2 = (th_1 + wave_max)/2; %high_thresh
        %Loop 'discretizes' values and ensures transitions aren't falsly made into
        %a descretized point. 
        prev_state = 0;
        cur_state = 0;
        state_ctr = 10;
        interp_factor = 1;
        for i = 1:length(wave)
            if wave(i) > high_thresh
                cur_state = .06;
            elseif (wave(i)<high_thresh && wave(i)>midpt)
                cur_state = .04;
            elseif (wave(i)<midpt && wave(i)>low_thresh)
                cur_state = .02;
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

        time_cutoff = 0.00; %Trims min bit periods found to remove errors
        time_differences = time_differences(time_differences>time_cutoff*3*41e-12); %Can cut out super low unexpected results
        time_differences = sort(time_differences);
        time_differences_norm = time_differences/min(time_differences(floor(0.5*length(time_differences)):end));

        bin = time_differences(time_differences_norm > 0.5 & time_differences_norm < 1.5);
        freq_recov = 1/mean(bin)
%         histogram(time_differences/1e-9, 500);
%         title("Distribution of Time Differences");
%         xlabel("Time Differences [ns]");
%         ylabel("Number of Occurences");
        
        freq_recov_save(p) = freq_recov;
        
        %         % Find waveform crossing points
%         pct20 = 0.2*eye_amp + zero_lvl_mean;
%         pct50 = (one_lvl_mean+zero_lvl_mean)/2;
%         pct80 = 0.8*eye_amp + zero_lvl_mean;
%         [crossing middlepts] = findcross(time_eye, eye, packet, [pct20, pct50, pct80]);
%         L = length(crossing);
% 
%         crossing_pct = (pct50 - zero_lvl_mean)/(one_lvl_mean - zero_lvl_mean);
% 
%         [crossing middlepts] = findcross(time_eye, eye, packet, [pct20, pct50, pct80]);
%         L = length(crossing);
% 
%         % First Crossing points
%         size(crossing(1, 1:2, 2, :));
%         crossT = reshape(crossing(1, 1:2, 2, :), [1, 2*L]);
% 
%         % Calculate CrossT mean and std
%         crossT_mean = mean(crossT);
%         crossT_std = std(crossT);
% 
%         jitter_rms = crossT_std;
%         jitter_pp = max(crossT)-min(crossT);
% 
%         % Second Crossing points
%         crossT2 = reshape(crossing(2, 1:2, 2, :), [1, 2*L]);
% 
%         % Calculate CrossT2 mean and std
%         crossT2_mean = mean(crossT2);
%         crossT2_std = std(crossT2);
% 
%         % Eye Width
%         eye_width = (crossT2_mean - 3*crossT2_std) - (crossT_mean + 3*crossT_std);
% 
%         % Rise time
%         rise_pct20_cross = reshape(crossing(1, 1, 1, :), [1, L]);
%         rise_pct80_cross = reshape(crossing(1, 1, 3, :), [1, L]);
% 
%         rise_time = rise_pct80_cross - rise_pct20_cross;
%         rise_time_mean = mean(rise_time);
%         rise_time_std = std(rise_time);
% 
%         % Fall time
%         fall_pct20_cross = reshape(crossing(1, 2, 1, :), [1, L]);
%         fall_pct80_cross = reshape(crossing(1, 2, 3, :), [1, L]);
%         fall_time = fall_pct20_cross - fall_pct80_cross;
%         fall_time_mean = mean(fall_time);
%         fall_time_std = std(fall_time);
    end

    % Plots
    if(plot_data == 1)
        h = figure();
        set(h,'WindowStyle','docked');

        if(wave_type == 3)
            surf(X, Y, scaled_eye');
            view(0,90);
            shading interp;
            xlim([time_eye_ps(1), time_eye_ps(end)]);
            cmap = colormap(hot);
%             cmap(1, :) = [159, 182, 205] ./256;
%             colormap(cmap);
            set(gca,'Color',cmap(1, :));
            set(gca,'visible','off');
            grid off;
            title("Unmodulated 2 GHz Signal and 4 GHz BW");
            xlabel("Time [ps]");
            ylabel("Voltage [mV]");
%             print -depsc high_eye
        else
            surf(X, Y, scaled_eye');
            view(0,90);
            shading interp;
            xlim([time_eye_ps(1), time_eye_ps(end)]);
            cmap = colormap(turbo);
%             colorbar;
            set(gca,'Color',cmap(1, :));
            title("0 km Fiber - Eye Diagram");
            xlabel("Time [ps]");
            ylabel("Voltage [mV]");
%             ylim([-1, 45])
%             print -depsc 0km
        end
    elseif(plot_data == 2)
%         % Plot time waveform
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plot(time, wave); hold on;
%         plot(time_interp, wave_interp);
%         title("Eye Diagram - 6GHz signal - 4 GHz Bandwidth");
%         xlabel("Time [ps]");
%         ylabel("Voltage [V]");
% 
%         % Plot eye diagram
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plot(time_eye_ps, eye_mv);
%         xlim([time_eye_ps(1), time_eye_ps(end)]);
%         title("Eye Diagram - 6GHz signal - 4 GHz Bandwidth");
%         xlabel("Time [ps]");
%         ylabel("Voltage [mV]");
%         yline(high_thresh/mV, 'linewidth', 2);
%         yline(low_thresh/mV, 'linewidth', 2);
%         xline(-ten_pct/ps);
%         xline(ten_pct/ps);
%         yline(midpt/mV, 'linewidth', 2);
% 
%         % Plot eye diagram gradient
%         h = figure();
%         set(h,'WindowStyle','docked');
%         surf(X, Y,scaled_eye');
%         view(0,90);
%         shading interp;
%         xlim([time_eye_ps(1), time_eye_ps(end)]);
%         cmap = colormap(turbo);
%         colorbar;
%         set(gca,'Color',cmap(1, :))
% 
%         % Plot one-level voltage distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(one_lvl_avg, 10)
%         title("One-Level Voltage Distribution")
%         xlabel("One-Level Voltage [V]")
%         ylabel("PDF")
%         legend("Experimental Data", "Gaussian Curve Fit")
% 
%         % Plot zero-level voltage distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(zero_lvl_avg, 10)
%         title("Zero-Level Voltage Distribution")
%         xlabel("Zero-Level Voltage Level [V]")
%         ylabel("PDF")
%         legend("Experimental Data", "Gaussian Curve Fit")
% 
%         % Plot zero-level voltage distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(mid0_lvl_avg, 10)
%         title("mid0 Voltage Distribution")
%         xlabel("mid0 Voltage Level [V]")
%         ylabel("PDF")
%         legend("Experimental Data", "Gaussian Curve Fit")
%         
%         % Plot zero-level voltage distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(mid1_lvl_avg, 10)
%         title("mid1 Voltage Distribution")
%         xlabel("mid1 Voltage Level [V]")
%         ylabel("PDF")
%         legend("Experimental Data", "Gaussian Curve Fit")
%         % Plot Jitter distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(crossT, 10)
%         title("Jitter Distribution");
%         xlabel("Time of 50% Crossing [s]");
%         ylabel("PDF");
%         legend("Experimental Data", "Gaussian Curve Fit")
% 
%         % Plot Rise Time distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(rise_time, 10)
%         title("Rise Time Distribution");
%         xlabel("Rise Time [ps]");
%         ylabel("PDF");
% 
%         % Plot Fall Time distribution
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plothist(fall_time, 10)
%         title("Fall Time Distribution");
%         xlabel("Fall Time [ps]");
%         ylabel("PDF");

%         h = figure();
%         set(h,'WindowStyle','docked');
%         endNbr = 10000;
%         plot(time_interp(1:endNbr)/ps, wave_interp(1:endNbr)/mV);
%         hold on;
%         plot(time_interp(65:endNbr)/ps, driveSig(1:endNbr-64)/mV);
%         plot(time_interp(65:endNbr)/ps, prbsSig(1:endNbr-64)/mV);
%         legend('original','recovered','prbs');
%         xlim([time_interp(1000)/ps, time_interp(2000)/ps]);
%         
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plot(bit_corr);
%         
%         h = figure();
%         set(h,'WindowStyle','docked');
%         plot(bit_stream); 
%         hold on;
%         plot(shifted_prbs);
%         xlim([1, 20]);
%         legend("recovered", "prbs");
%         
%         h = figure();
%         set(h,'WindowStyle','docked');     
%         plot(bit_stream_prbs);
%         hold on;
%         plot(sumSig);
%         xlim([1, 40]);
%         legend("bits", "sumSig");
    elseif(plot_data == 3)
        if(wave_type < 4)
        
            subplot(3, 3, p)
            surf(X, Y, scaled_eye');
            view(0,90);
            shading interp;
            xlim([time_eye_ps(1), time_eye_ps(end)]);
            if(wave_type == 1)
                ylim([-5, 45]);
            end
            cmap = colormap(turbo);
            colorbar;
            set(gca,'Color',cmap(1, :));
            if(wave_type == 0)
                plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
            elseif(wave_type == 1)
                plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
            elseif(wave_type == 2)
                plot_title = strcat(num2str((floor((p-1)/3)+1)*3), "km Fiber + ", num2str(2^(2-mod(p-1, 3))), "GHz BW");
            else
                plot_title = "DEFAULT";
            end
            title(plot_title);
            xlabel("Time [ps]");
            ylabel("Voltage [mV]");

            if(wave_type == 1)
                sgtitle('Eye Diagrams from varying Modulated Signal Frequencies');
            elseif(wave_type == 2)
                sgtitle('Eye Diagrams from varying Fiber Lengths');
            end
        else
%             row = mod(p-1, 3);
%             col = ; 
%             pos = row*7 + col
            index = mod(p-1, 3)*6 + mod(p-1, 7) + 1
            
            if(mod(p-1, 3) == 0)
                subplot(3, 7, a)
                a = a+1;
                surf(X, Y, scaled_eye');
                view(0,90);
                shading interp;
                xlim([time_eye_ps(1), time_eye_ps(end)]);
                if(wave_type == 1)
                    ylim([-5, 45]);
                end
                cmap = colormap(turbo);
                colorbar;
                set(gca,'Color',cmap(1, :));
                if(wave_type == 0)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 1)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 2)
                    plot_title = strcat(num2str((floor((p-1)/3)+1)*3), "km Fiber + ", num2str(2^(2-mod(p-1, 3))), "GHz BW");
                else
                    plot_title = p;
                end
                title(plot_title);
                xlabel("Time [ps]");
                ylabel("Voltage [mV]");

                if(wave_type == 1)
                    sgtitle('Eye Diagrams from varying Modulated Signal Frequencies');
                elseif(wave_type == 2)
                    sgtitle('Eye Diagrams from varying Fiber Lengths');
                end
            elseif(mod(p-1, 3) == 1)
                subplot(3, 7, b+7)
                b = b+1;
                surf(X, Y, scaled_eye');
                view(0,90);
                shading interp;
                xlim([time_eye_ps(1), time_eye_ps(end)]);
                if(wave_type == 1)
                    ylim([-5, 45]);
                end
                cmap = colormap(turbo);
                colorbar;
                set(gca,'Color',cmap(1, :));
                if(wave_type == 0)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 1)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 2)
                    plot_title = strcat(num2str((floor((p-1)/3)+1)*3), "km Fiber + ", num2str(2^(2-mod(p-1, 3))), "GHz BW");
                else
                    plot_title = p;
                end
                title(plot_title);
                xlabel("Time [ps]");
                ylabel("Voltage [mV]");

                if(wave_type == 1)
                    sgtitle('Eye Diagrams from varying Modulated Signal Frequencies');
                elseif(wave_type == 2)
                    sgtitle('Eye Diagrams from varying Fiber Lengths');
                end
            elseif(mod(p-1, 3) == 2)
                subplot(3, 7, c+14)
                c = c+1;
                surf(X, Y, scaled_eye');
                view(0,90);
                shading interp;
                xlim([time_eye_ps(1), time_eye_ps(end)]);
                if(wave_type == 1)
                    ylim([-5, 45]);
                end
                cmap = colormap(turbo);
                colorbar;
                set(gca,'Color',cmap(1, :));
                if(wave_type == 0)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 1)
                    plot_title = strcat(num2str(sig_freq), "GHz Signal + ", num2str(4-mod(p-1, 3)), "GHz BW");
                elseif(wave_type == 2)
                    plot_title = strcat(num2str((floor((p-1)/3)+1)*3), "km Fiber + ", num2str(2^(2-mod(p-1, 3))), "GHz BW");
                else
                    plot_title = p;
                end
                title(plot_title);
                xlabel("Time [ps]");
                ylabel("Voltage [mV]");

                if(wave_type == 1)
                    sgtitle('Eye Diagrams from varying Modulated Signal Frequencies');
                elseif(wave_type == 2)
                    sgtitle('Eye Diagrams from varying Fiber Lengths');
                end
            end
        end
    end
end

%% Plot BER vs Connected Fiber Length
h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
semilogy(fiber_L, BER_save(1:3:end), '-*', 'linewidth', 2);
hold on;
semilogy(fiber_L, BER_save(2:3:end), '-*', 'linewidth', 2);
semilogy(fiber_L, BER_save(3:3:end), '-*', 'linewidth', 2);
legend("0.975 GHz", "1.95 GHz", "3.90 GHz", 'location', 'southeast');
xticks(fiber_L);
xlabel("Fiber Length [km]");
ylabel("BER");
title({"Recovered BER vs. Connected Fiber Length","at Varying Frequencies"});
% print -depsc BER_vs_fiberLength
% exportgraphics(h,'BER_vs_fiberLength.png','Resolution',300)

%% Plot variance over fiber length
close all; clc;

frq = 3;

h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
plot(fiber_L, one_lvl_mean_norm_save(frq:3:end), '-*', 'linewidth', 2);
hold on;
plot(fiber_L, mid1_lvl_mean_norm_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, mid0_lvl_mean_norm_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, zero_lvl_mean_norm_save(frq:3:end), '-*', 'linewidth', 2);
legend("1", "2/3","1/3", '0', 'location', 'southeast');
xticks(fiber_L);

h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
plot(fiber_L, one_lvl_std_norm_save(frq:3:end), '-*', 'linewidth', 2);
hold on;
plot(fiber_L, mid1_lvl_std_norm_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, mid0_lvl_std_norm_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, zero_lvl_std_norm_save(frq:3:end), '-*', 'linewidth', 2);
legend("1", "2/3","1/3", '0', 'location', 'southeast');
xticks(fiber_L);

% h = figure();
% set(h,'WindowStyle','docked');
% fiber_L = [0, 3, 6, 9, 12, 15, 25];
% plot(fiber_L, one_lvl_std_norm_save(1:3:end), '-*', 'linewidth', 2);
% hold on;
% plot(fiber_L, one_lvl_std_norm_save(2:3:end), '-*', 'linewidth', 2);
% plot(fiber_L, one_lvl_std_norm_save(3:3:end), '-*', 'linewidth', 2);
% plot(fiber_L, zero_lvl_std_norm_save(1:3:end), '-*', 'linewidth', 2);
% plot(fiber_L, zero_lvl_std_norm_save(2:3:end), '-*', 'linewidth', 2);
% plot(fiber_L, zero_lvl_std_norm_save(3:3:end), '-*', 'linewidth', 2);
% legend("0.975 GHz", "1.95 GHz", "3.90 GHz", 'location', 'southeast');
% xticks(fiber_L);

h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
plot(fiber_L, Q_01_save(frq:3:end), '-*', 'linewidth', 2); hold on;
plot(fiber_L, Q_02_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, Q_03_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, Q_12_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, Q_13_save(frq:3:end), '-*', 'linewidth', 2);
plot(fiber_L, Q_23_save(frq:3:end), '-*', 'linewidth', 2);
legend("Q_{01}","Q_{02}","Q_{03}","Q_{12}","Q_{13}","Q_{23}", 'location', 'northeast');
xticks(fiber_L);
title("Q Value for Each Level Eye over Fiber Length");
% print -depsc Q_vs_fiberLength
% exportgraphics(h,'Q_vs_fiberLength.png','Resolution',300);

%% Freq Recovery
frq = 1;
GHz = 1e9;
h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
plot(fiber_L, freq_recov_save(1:3:end)/GHz, '-*', 'linewidth', 2);
hold on;
plot(fiber_L, freq_recov_save(2:3:end)/GHz, '-*', 'linewidth', 2);
plot(fiber_L, freq_recov_save(3:3:end)/GHz, '-*', 'linewidth', 2);
yline(0.975);
yline(1.95);
yline(3.9);
legend("0.975 GHz", "1.95 GHz", "3.90 GHz", 'location', 'southeast');
xticks(fiber_L);

%% Insertion Loss
h = figure();
set(h,'WindowStyle','docked');
fiber_L = [0, 3, 6, 9, 12, 15, 25];
p = fit(fiber_L', one_lvl_mean_save(1:3:end)', 'exp1')
% y = polyval(p, fiber_L);

semilogy(fiber_L, one_lvl_mean_save(1:3:end)/mV, '-*', 'linewidth', 2);hold on;
semilogy(fiber_L, one_lvl_mean_save(2:3:end)/mV, '-*', 'linewidth', 2);
semilogy(fiber_L, one_lvl_mean_save(3:3:end)/mV, '-*', 'linewidth', 2);
semilogy(fiber_L, p.a*exp(p.b*fiber_L)/mV, 'linewidth', 2);
legend("0.975 GHz", "1.95 GHz", "3.90 GHz", 'Fitted Line', 'location', 'northeast');
xticks(fiber_L);
xlabel("Fiber Length [km]");
ylabel("One Level Mean Voltage [mV]");
title("Voltage vs. Fiber Length");

print -depsc V_vs_fiberLength
exportgraphics(h,'V_vs_fiberLength.png','Resolution',300);
