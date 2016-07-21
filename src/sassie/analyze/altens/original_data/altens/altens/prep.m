function [Di,vcoor,rlist]=prep(D,NH,reslist)
%-------------------------------------------------------
%  df-oct-98 (modified/corrected version) df-may-98
%	prepare working set of data 
%	(for selected residues)
%	to be used by r2r1fit
%
%	uses: combine2.m
%------------------------------------------------------- 

nNH=size(NH,1);
%--------make res-list--------
if isempty(reslist), 
   reslist=sort(NH(:,1)); 
else
   reslist=sort(reslist(:)); 
end
nres0=length(reslist);
%----------- prepare & normalize vNH ------------
vNH=NaN*ones(nres0,4);
vNH(:,1)=reslist;
for ii=1:nres0,
    ind=find(NH(:,1)==vNH(ii,1));
    if ~isempty(ind),
        vNH(ii,2:4)=NH(ind,2:4)/sqrt(NH(ind,2:4)*NH(ind,2:4)');
    end
end
%---------- correct for HF-components----------

 vDi=NaN*ones(nres0,2);
 vDi(:,1)=reslist;
 for ii=1:nres0,
    ind=find(D(:,1)==vDi(ii,1));
    if ~isempty(ind),
        vDi(ii,2:3)=D(ind,2:3);
    end
end

%--------- select residues-------
sel=find((~isnan(vNH(:,2))).*(~isnan(vDi(:,2))));
%----------select working set of data-----------
Di=vDi(sel,2:3);
vcoor=vNH(sel,2:4);
rlist=vNH(sel,1);
return
%===============================================