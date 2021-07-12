
function [S,A,Do]=RDCcalc(RDC,rNH)
%--------------------------------------------
%   convert principal value of the RDC tensor
%   into alignment (A) and Saupe (S) tensor
% ------------------------------------------- 
if nargin < 2 rNH=1.04; end
gH = 26750.6;
gN = -2710.6;
rNH = rNH*1e-8;
h = 6.6262e-27;
Do = -gH*gN*(h/2/pi)/pi/rNH^3;      %strength of DipCoupling
Do = 2.17e4;
S = RDC/Do;                         %Saupe tensor
A = 2/3*S;                           %Aligment tensor

%calculating static dipolar constant
% Do=(uo/2pi)gH*gN*(h/2pi)/pi/rNH^3
%in SGS (uo/2pi)=1
