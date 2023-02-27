function [ratio,sigma,vcoor,rlist]=r2r1prep(r2,r1,NH,reslist,r3)
%-------------------------------------------------------
%  df-feb-12 modified error estimation for ratio
%  df-oct-98 (modified/corrected version) df-may-98
%	prepare working set of data 
%	(for selected residues)
%	to be used by r2r1fit
%
%	uses: combine2.m
%------------------------------------------------------- 
fR2=(0.87/0.955)^2;				%modified
fR1=(0.87/0.921)^2;
fRj=(0.87/1)^2;
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
if ~isempty(r3),
  g2g=2710.5/26750;

  z13=combine2(r1,r3,reslist,2);
  z13er=combine2(r1,r3,reslist,3);
  Jh=[z13(:,1),z13(:,2).*(1-z13(:,3))*g2g/5];
  Jh=[Jh,Jh(:,2).*sqrt((z13er(:,2)./z13(:,2)).^2+(z13er(:,3)./(1-z13(:,3))).^2)];
  Jher=Jh(:,3);
  
  NOE=z13(:,3);             %df-12
  NOEer=z13er(:,3);         %df-12
  
  z1=combine2(r1,Jh,reslist,2);
  z1er=combine2(r1,Jh,reslist,3);
  r1a=z1(:,2)-z1(:,3)*7*fR1;
  r1er=z1er(:,2);
  
  R1=z1(:,2);               %df-12
  
  z2=combine2(r2,Jh,reslist,2);
  z2er=combine2(r2,Jh,reslist,3);
  r2a=z2(:,2)-z2(:,3)*13/2*fR2;
  r2er=z2er(:,2);

  R2=z2(:,2);               %df-12
 
  %-------------calc. RJ=2R2adj-R1adj--------------
  zj=combine2(r2,r1,reslist,2);			
  zjer=combine2(r2,r1,reslist,3);
  yj=combine2(r2,Jh,reslist,2);
  yjer=combine2(r2,Jh,reslist,3);
  rj=2*zj(:,2)-zj(:,3)-6*yj(:,3)*fRj;   %modified
  
  Rj=rj;                    %df-12
  
  C1=7/5*g2g*fR1;           %df-12
  C3=6/5*g2g*fRj;           %df-12
  
  %--------------calc. RATIO=R1adj/RJ--------------
  RATIO=r1a./rj;
  %-------------calc. RATIOerror -- corrected------
  RATIOer=(R1./Rj).*(R2./Rj).*sqrt(4*((r1er./R1).^2+(r2er./R2).^2).*...
      (1-C1*(1-NOE))+(NOEer.^2).*(2*C1-(R1./R2).*(C1+C3)).^2);          %df-12
   
%  RATIOer1=RATIO.*sqrt((r1er.*(1./r1a+1./rj)).^2+...
%     (2*r2er./rj).^2+(Jher.*(7*fR1./r1a-6*fRj./rj)).^2);
% [RATIO RATIOer RATIOer1]
else				%NOE data are not available
  z1=combine2(r2,r1,reslist,2);
  z1er=combine2(r2,r1,reslist,3);
  rj=2*z1(:,2)-z1(:,3);
  r1=z1(:,3);
  r1er=z1er(:,3); 
  r2er=z1er(:,2); 
  RATIO=r1./rj;
  RATIOer=RATIO.*sqrt((r1er.*(1./r1+1./rj)).^2+...
     (2*r2er./rj).^2);
end 

%--------- select residues-------
sel=find((~isnan(vNH(:,2))).*(~isnan(RATIO)));
%----------select working set of data-----------
ratio=RATIO(sel);
sigma=RATIOer(sel);
vcoor=vNH(sel,2:4);
rlist=z1(sel,1);
return
%===============================================