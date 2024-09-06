function f=calc_Dz(r,Tx);


%calculates Dzz from ratio and TAUx

f=(r/(6*Tx*1e-9))*1e-7;
return