clearvars; clc;
close all;

% SETUP
% 0 = No plots, Calculate param 
% 1 = Plot Gradient Eye Only, No param calculations, Used to center graph
% 2 = Plot everything, Calculate all param
% 3 = Plot all traces in Grid
plot_data = 1;

% -1-select fiber
% 0-no fiber, no mod
% 1-no fiber, mod
% 2-fiber, mod
% 3-lab4 day1, wrong vPI
wave_type = -1;
select = 1;

% Initialize Data Arrays
Q_save = zeros(1, 9);
BER_save = zeros(1, 9);

if(select)
    nbrTraces = 1:1;
elseif(wave_type > 0)
    nbrTraces = 1:9;
else
    nbrTraces = 1:8;
end

if(plot_data == 3)
    h = figure();
    h.Position = [100 100 1500 900];
%     set(h,'WindowStyle','docked');
end

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
        filename = strcat('media_day2/lab4_002',selectdata,'_ALL.csv');
    else
        filename = '\media\lab5_day4_004_ALL.csv';
    end
    T = readtable(filename);
    data = table2array(T);

    % Separate data to time and wave vectors
    time = data(6:end, 1);
    if(wave_type > 0)
        wave = data(6:end, 3);
    else 
        wave = data(6:end, 2);
    end

    % 2 GHz = .5 ns  - 500ps, 4000 pulses in 2us
    % 50e3pts, 40ps/pt, 12.5 pts per pulse,
    % interp - 400e3 pts, 20ps/pt, 25pts per pulse
    if(wave_type == -1)
        sig_freq = 1.95;
    elseif(wave_type < 2)
        sig_freqs = [2, 2, 2, 3, 3, 3, 4, 4, 4];
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
    else
        sig_offsets = 5;
    end
    
    offset = sig_offsets(p);
    
    % Waveform properties
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

        % Calculate I0 mean and std
        zero_lvl_mean = mean(zero_lvl_avg);
        zero_lvl_std = std(zero_lvl_avg);
 
        mid0_lvl_mean = mean(mid0_lvl_avg);
        mid0_lvl_std = std(mid0_lvl_avg);
        
        mid1_lvl_mean = mean(mid1_lvl_avg);
        mid1_lvl_std = std(mid1_lvl_avg); 

        % Calculate Eye Amplitude
        eye_amp = one_lvl_mean - zero_lvl_mean;
        eye_ht = one_lvl_mean - 3*one_lvl_std - zero_lvl_mean - 3*zero_lvl_std;

        % Calculate Q value and BER
        Q = eye_amp/(one_lvl_std+zero_lvl_std);
        BER = 1/2 * erfc(Q./sqrt(2));
        
        Q_save(p) = Q';
        BER_save(p) = BER';

        %Find New wave crossing points
        midpt = (mid0_lvl_std*mid1_lvl_mean + mid1_lvl_std*mid0_lvl_mean)...
            /(mid0_lvl_std+mid1_lvl_std);
        high_thresh = (mid1_lvl_std*one_lvl_mean + one_lvl_std*mid1_lvl_mean)...
            /(one_lvl_std+mid1_lvl_std);
        low_thresh = (mid0_lvl_std*zero_lvl_mean + zero_lvl_std*mid0_lvl_mean)...
            /(zero_lvl_std+mid0_lvl_std);
        
        % Recover Bit Stream
        avg_middle = mean(middledata);
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
        
        O = 15;
        offset = 6;

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
        
        find((shifted_prbs - bit_stream) ~= 0)
        BER = length(find((shifted_prbs - bit_stream) ~= 0))/length(bit_stream)
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

        h = figure();
        set(h,'WindowStyle','docked');
        endNbr = 10000;
        plot(time_interp(1:endNbr)/ps, wave_interp(1:endNbr)/mV);
        hold on;
        plot(time_interp(65:endNbr)/ps, driveSig(1:endNbr-64)/mV);
        plot(time_interp(65:endNbr)/ps, prbsSig(1:endNbr-64)/mV);
        legend('original','recovered','prbs');
        xlim([time_interp(1000)/ps, time_interp(2000)/ps]);
        
        h = figure();
        set(h,'WindowStyle','docked');
        
        plot(bit_corr);
        
        h = figure();
        set(h,'WindowStyle','docked');
        
        plot(bit_stream); 
        hold on;
        plot(shifted_prbs);
        xlim([1, 20]);
        legend("recovered", "prbs")
        
    elseif(plot_data == 3)
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
    end
end
