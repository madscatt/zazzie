% computing dipolar constants
mu0 = 4.0*pi*1.0e-7;
h = 6.62606896e-34;
gH = 42.576e6*2*pi;   %Hz/T  or 267.513 *1e6 rad/s/T
gN = -4.3156e6*2*pi;  %Hz/T  or -27.116 *1e6 rad/s/T
gC = 10.705e6*2*pi;   %Hz/T  or 67.262  *1e6 rad/s/T

rNH=1.04e-10;
rCH=1.09e-10;

C_nh = mu0*abs(gN)*gH*h/(2*pi)^3/(rNH)^3
C_ch = mu0*gC*gH*h/(2*pi)^3/(rCH)^3


%PATI: C_nh = 1.836e+04, C_ch = 3.731e+04