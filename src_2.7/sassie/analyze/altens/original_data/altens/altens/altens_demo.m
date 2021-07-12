%altens_demo -- demo script for running altens
%it uses not the best data set, but for demo purposes it should
%suffice

load demo
nsteps = 500;
[table_all,Ds]=altens(reslist,vNH,D,nsteps);
