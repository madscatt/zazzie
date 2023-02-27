function z=r2r1mc_iso(par,wN)
tauc=par;
wT2=wN*wN*tauc*tauc;


z=-(3/4)*((2*tauc*wN^2)/(1+wT2)^2);

z=z'*z;
return