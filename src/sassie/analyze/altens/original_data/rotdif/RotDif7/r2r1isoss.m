function z=r2r1isoss(param,ratio,sigma,wN)
%-------------------------------------------
%  df-oct-98 an economical version   df-may-98
%	calc.ss of (R2R1expt-R2R2cal)/sigma)
%	param=[Tx,Dz/Dx]
%-------------------------------------------
Tx=param(1);
wT2=wN*wN*Tx*Tx;
%------------calculate ratio ---------
calc=3/4/(1+wT2);
%------------get difference-----------
diff=(ratio-calc)./sigma;
z=diff'*diff;				%might need looping for MCC
return
%===============================================
