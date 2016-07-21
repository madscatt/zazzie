function Q=quality_iso(param,ratio,sigma,wN)
%-------------------------------------------
%  df-oct-98 an economical version   df-may-98
%	calc.ss of (R2R1expt-R2R2cal)/sigma)
%	param=[Tx,Dz/Dx]
%-------------------------------------------
Tx=param(1);
wT2=wN*wN*Tx*Tx;
%------------calculate ratio ---------
calc=3/4/(1+wT2);
lr=size(ratio,1);
calc=ones(lr,1).*calc;
%------------quality factor-----------

num=(ratio-calc).^2;
num2=mean(num);
den=ratio-mean(ratio);
den1=den.^2;
den2=mean(den1);
den3=2*den2;
Q=sqrt(num2/den3);





return
%===============================================
