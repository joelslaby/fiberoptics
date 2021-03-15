function [crossing, middlepts] = findcross(time_eye, eye, packet, thresh)
    %time_eye, eye, first, thresh, rise, packet
    % Dimensions [1st/2nd, Rise/Fall, Thresh, Crossings]
    order = 0; %1st=1, 2nd=2
    traj = 0; % Rise=1, Fall=2
    pct = 0; % 20%=1 50%=2 80%=3
    L = 2*length(eye);
    crossing = zeros(2, 2, 3, L);
    ctr = ones(2, 2);
    
    minrise = [0, 0, 0; 0, 0, 0];
    maxrise = [0, min(min(eye)), 0; 0, min(min(eye)), 0];
    minfall = [0, 0, 0; 0, 0, 0];
    maxfall = [0, min(min(eye)), 0; 0, min(min(eye)), 0];
%     [L, ~] = size(eye);
    
    for i=1:length(eye)
        wave = eye(:, i)';
        if(wave(1)<thresh(2)) % Checks - to + Rise
            traj = 1;
            if((find(wave>thresh(2), 1) < packet)) % Checks 1st or 2nd
                order = 1;
            elseif((find(wave>thresh(2), 1) > packet))
                order = 2;
            else
                order = nan;
            end
            if(~isnan(order))
                for k=1:length(thresh)
                    loc = find(wave>thresh(k), 1);
                    if(loc == 1)
                        loc = [];
                    end
                    if(~isempty(loc))
                        cross = InterX([time_eye(loc), time_eye(loc-1); wave(loc), wave(loc-1)], [time_eye(loc), time_eye(loc-1); thresh(k), thresh(k)]);
                        if(cross(1) < minrise(order, 2) & k == 2)
                            minrise(order, 1) = i;
                            minrise(order, 2) = cross(1);
                            minrise(order, 3) = cross(2);
                        elseif(cross(1) > maxrise(order, 2) & k == 2)
                            maxrise(order, 1) = i;
                            maxrise(order, 2) = cross(1);
                            maxrise(order, 3) = cross(2);
                        end
                        crossing(order, traj, k, ctr(traj, order)) = cross(1);
                    end
                end
                ctr(traj, order) = ctr(traj, order) + 1;
            end 
        elseif(wave(1)>thresh(2)) % Checks + to - Fall
            traj = 2;
            if((find(wave<thresh(2), 1) < packet)) % Checks 1st or 2nd
                order = 1;
            elseif((find(wave<thresh(2), 1) > packet))
                order = 2;
            else
                order = nan;
            end
            if(~isnan(order))
                for k=1:length(thresh)
                    loc = find(wave<thresh(k), 1);
                    if(loc == 1)
                        loc = [];
                    end
                    if(~isempty(loc))
                        cross = InterX([time_eye(loc), time_eye(loc-1); wave(loc), wave(loc-1)], [time_eye(loc), time_eye(loc-1); thresh(k), thresh(k)]);
                        if(cross(1) < minfall(order, 2) & k == 2)
                            minfall(order, 1) = i;
                            minfall(order, 2) = cross(1);
                            minfall(order, 3) = cross(2);
                        elseif(cross(1) > maxfall(order, 2) & k == 2)
                            maxfall(order, 1) = i;
                            maxfall(order, 2) = cross(1);
                            maxfall(order, 3) = cross(2);
                        end
                        crossing(order, traj, k, ctr(traj, order)) = cross(1);
                    end
                end
                ctr(traj, order) = ctr(traj, order) + 1;
            end 
        end
    end
    [minrise(1, 2), maxrise(1, 2), minfall(1, 2), maxfall(1, 2)];
    middlepts = [minrise(1,:); maxrise(1,:); minfall(1,:); maxfall(1,:)];
%     for i=1:2
%         for j=1:2
%             for k=1:3
%                 c = ctr(i, j)
%                 crossing(i, j, k, c:L) = nan;
%             end
%         end
%     end
    crossing = crossing(1:2, 1:2, 1:3, 1:ctr-1);
%     crossing = crossing(~isnan(crossing));

%     crossing = [];
%     if(first)
%         for i=1:length(eye)
%             if((eye(1, i)<thresh) & (find(eye(:, i)>thresh, 1) < packet) & rise~=2)
%                 id1 = find(eye(:, i)>thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];
%             elseif((eye(1, i)>thresh) & (find(eye(:, i)<thresh, 1) < packet) & rise~=1)
%                 id1 = find(eye(:, i)<thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];   
%             end
%         end
%     else
%         for i=1:length(eye)
%             if((eye(1, i)<thresh) & (find(eye(:, i)>thresh, 1) > packet) & rise~=2 )
%                 id1 = find(eye(:, i)>thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];
%             elseif((eye(1, i)>thresh) & (find(eye(:, i)<thresh, 1) > packet) & rise~=1)
%                 id1 = find(eye(:, i)<thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];   
%             end
%         end
%     end
end