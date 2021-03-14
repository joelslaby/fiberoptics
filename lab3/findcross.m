function [crossing] = findcross(time_eye, eye, first, thresh, rise)
    crossing = [];
%     if(rise)
%         for i=1:length(eye)
%             if((eye(1, i)<0) & (find(eye(:, i)>0, 1) < 25))
%                 id1 = find(eye(:, i)>thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];
%             end
%         end
%     else
%         for i=1:length(eye)
%             if((eye(1, i)>0) & (find(eye(:, i)<0, 1) < 25))
%                 id1 = find(eye(:, i)<thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];  
%             end
%     end
% 
%     if(first & ~rise)
%         crossing = [];
%         for i=1:length(eye)
%             if((eye(1, i)<0) & (find(eye(:, i)>0, 1) < 25))
%                 id1 = find(eye(:, i)>thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];
%             elseif((eye(1, i)>0) & (find(eye(:, i)<0, 1) < 25))
%                 id1 = find(eye(:, i)<thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];   
%             end
%         end
%     else
%         crossing = [];
%         for i=1:length(eye)
%             if((eye(1, i)<0) & (find(eye(:, i)>thresh, 1) > 25))
%                 id1 = find(eye(:, i)>thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];
%             elseif((eye(1, i)>0) & (find(eye(:, i)<thresh, 1) > 25))
%                 id1 = find(eye(:, i)<thresh, 1);
%                 cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
%                 crossing = [crossing, cross(1)];   
%             end
%         end
%     end
    
    if(first)
        crossing = [];
        for i=1:length(eye)
            if((eye(1, i)<0) & (find(eye(:, i)>0, 1) < 25) & rise~=2)
                id1 = find(eye(:, i)>thresh, 1);
                cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
                crossing = [crossing, cross(1)];
            elseif((eye(1, i)>0) & (find(eye(:, i)<0, 1) < 25) & rise~=1)
                id1 = find(eye(:, i)<thresh, 1);
                cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
                crossing = [crossing, cross(1)];   
            end
        end
    else
        crossing = [];
        for i=1:length(eye)
            if((eye(1, i)<0) & (find(eye(:, i)>thresh, 1) > 25) & rise~=2 )
                id1 = find(eye(:, i)>thresh, 1);
                cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
                crossing = [crossing, cross(1)];
            elseif((eye(1, i)>0) & (find(eye(:, i)<thresh, 1) > 25) & rise~=1)
                id1 = find(eye(:, i)<thresh, 1);
                cross = InterX([time_eye(id1), time_eye(id1-1); eye(id1, i), eye(id1-1, i)], [time_eye(id1), time_eye(id1-1); thresh, thresh]);
                crossing = [crossing, cross(1)];   
            end
        end
    end
     
end