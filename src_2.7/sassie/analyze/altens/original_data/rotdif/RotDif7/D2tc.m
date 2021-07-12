function g=D2tc(v1,v2,v3)


% calculates TAUc from Dxx, Dyy, Dzz

trace=v1+v2+v3;
g=1/(2*trace)*1e2;

