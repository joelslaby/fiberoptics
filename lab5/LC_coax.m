e0 = 8.854187e-12;
eR = 1;
u0 = 4*pi*1e-7;
uR = 1;

D = 0.81;
d = 5;

C = 2*pi*e0*eR/(log(D/d));
L = u0*uR*log(D/d)/(2*pi);

delay = sqrt(L*C)