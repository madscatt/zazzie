function x=sort_unused(reslist,rlist)
% extract unused residues
reslist=reslist;        %df modif
sl=length(rlist);       %df modif

tab=NaN*ones(sl,1);     %df modif
for ii=1:sl
  if (find(reslist==rlist(ii,1)))   %df modif
     tab(ii)=ii;        %df modif
  end
end

ind=find(isnan(tab));
x=rlist(ind);