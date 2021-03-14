clc; clearvars; close all;

%prelab 1 - p3

r1 = 4.1e-6;
r2 = 5.2e-6;
a1 = pi*r1^2;
a2 = pi*r2^2;
%a1 = r1;
%a2 = r2;
a3 = a2-a1;

Ng = 1.4682;
deln = .0036;

n1 = Ng * ((a1+a3*(1-deln))/(a1+a3))^(-1);
n2 = n1*(1-deln);

n2 = 1.4440;
n1 = n2/(1-deln)


NA = sqrt(n1^2-n2^2);
v = 2.405;

lambda = 2*pi/v * r1 * NA
