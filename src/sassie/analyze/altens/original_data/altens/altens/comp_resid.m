function [Di,new_ratio,rlist]=comp_resid(ratio,D)
%-------------------------------------------------------
%  ow-august 2002
%  compare ratio obtained by relxation data and res. dip. couplings
%  (for selected residues)
%  to be used 
%
%------------------------------------------------------- 


%----- make reslist ------------
if (size(ratio,1)>size(D,1))
    reslist=ratio(:,1);
end

if (size(D,1)>size(ratio,1))
    reslist=D(:,1);
end

if (size(D,1)==size(ratio,1))
    reslist=D(:,1);
end

reslist
nres0=length(reslist);

%----------- prepare ratio ------------
v_ratio=NaN*ones(nres0,4);
v_ratio(:,1)=reslist;
for ii=1:nres0,
    ind=find(ratio(:,1)==v_ratio(ii,1));
    if ~isempty(ind),
        v_ratio(ii,2)=ratio(ind,2);
    end
end
%---------- correct for HF-components----------

 vDi=NaN*ones(nres0,2);
 vDi(:,1)=reslist;
 for ii=1:nres0,
    ind=find(D(:,1)==vDi(ii,1));
    if ~isempty(ind),
        vDi(ii,2)=D(ind,2);
    end
end

%--------- select residues-------
sel=find((~isnan(v_ratio(:,2))).*(~isnan(vDi(:,2))));
%----------select working set of data-----------
Di=vDi(sel,2);
new_ratio=v_ratio(sel,2);
rlist=v_ratio(sel,1);




%------ plot ratio versus res. dip. couplings ----------
subplot(2,2,1)
plot(new_ratio,Di,'.r');
title('\fontsize{12}ratio versus res. dip. couplingfs');
xlabel('\fontsize{12}exp ratio');
ylabel('\fontsize{12}res. dip. couplings');




div=Di./new_ratio;
subplot(2,2,2)
plot(rlist,div,'.b');
title('\fontsize{12}residue number versus (D/ratio)');
xlabel('\fontsize{12}res. number');
ylabel('\fontsize{12}(D/ratio)');


div_abs=abs(Di)./new_ratio;
subplot(2,2,3)
plot(rlist,div_abs,'.g');
title('\fontsize{12}residue number versus (abs(D)/ratio)');
xlabel('\fontsize{12}res. number');
ylabel('\fontsize{12}(abs(D)/ratio)');


subplot(2,2,4)
plot(new_ratio,abs(Di),'.r');
title('\fontsize{12}ratio versus abs(D)');
xlabel('\fontsize{12}exp ratio');
ylabel('\fontsize{12}abs(D)');






return
%===============================================