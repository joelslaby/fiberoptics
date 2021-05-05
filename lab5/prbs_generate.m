close all; clearvars; clc;

O = 4;
offset = 6;

L = 2^(O)-1;

sig = prbs(O,L);
% sig = sig*2-1;
sigP = -0.5 * circshift(sig, -offset);
sumSig = (sig+sigP);

bit_stream_prbs = zeros(1, length(sig));
bit_stream_prbs(find(sumSig == -1.5)) = -1;
bit_stream_prbs(find(sumSig == -.5)) = -1/3;
bit_stream_prbs(find(sumSig == .5)) = 1/3;
bit_stream_prbs(find(sumSig == 1.5)) = 1;

N = 2^(O+1)-1;
temp = prbs(O,N);
corrMatrix = zeros(1, L+1);
summed = zeros(L, 1);
% for i = 1:L+1
%     sig2 = double(not(temp(1+i:L+i)));
%     R = corrcoef(sig1, sig2);
%     corrMatrix(i) = R(1, 2);
% end

% i = 1;
% sig2 = double(not(temp(1+i:L+i)))
% R = corrcoef(sig1, sig2)
% corrMatrix(i) = R(1, 2);

h = figure();
set(h,'WindowStyle','docked');

% off = 6;
% summed = sig + circshift(sig, off);

sig_short = sig(3:15)
sig
circshift(sig, -2);
out = conv(sig, flip(sig_short));
t = (1:length(out))-(length(out)+1);
plot(out);

% prbs(3, 2)

% h = figure();
% set(h,'WindowStyle','docked');
% 
% plot(t, conv(summed, flip(summed)));
