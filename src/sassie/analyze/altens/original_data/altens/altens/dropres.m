function newlist=dropres(rlist,EXCLUDE)
%---------------------------------------
%  df-oct-98
%	exclude residues from residue list
%---------------------------------------
nres0=length(rlist);
nexcl=length(EXCLUDE);
ind=ones(nres0,1);
for ii=1:nres0,
   if ~isempty(find(EXCLUDE(:)==rlist(ii))),
      ind(ii)=0;
   end
end
newlist=rlist(find(ind));
return
%=======================================
