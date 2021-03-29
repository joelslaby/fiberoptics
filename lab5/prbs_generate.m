close all; clearvars; clc;

O = 4;
off = 5;

L = 2^O-1;
N = 2^(O+1)-1;
sig1 = prbs(O,L);
temp = prbs(O,N);
corrMatrix = zeros(1, L+1);
summed = zeros(L, 1);
for i = 1:L+1
    sig2 = double(not(temp(1+i:L+i)));
    R = corrcoef(sig1, sig2);
    corrMatrix(i) = R(1, 2);
end

% i = 1;
% sig2 = double(not(temp(1+i:L+i)))
% R = corrcoef(sig1, sig2)
% corrMatrix(i) = R(1, 2);

