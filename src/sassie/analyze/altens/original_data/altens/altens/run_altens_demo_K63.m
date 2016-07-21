%altens_demo -- demo script for running altens
%it uses not the best data set, but for demo purposes it should
%suffice

load altens_demo_K63

save RDC.txt RDC_P_K63 -ascii
save vNH.txt vNH -ascii
save reslist.txt reslist -ascii


%[table_all,Ds,R,corr,tens,S_tensor]=altens(RDC_P_K63,vNH,reslist);
