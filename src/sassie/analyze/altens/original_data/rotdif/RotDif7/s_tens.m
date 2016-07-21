function SS=s_tens(vects,kfull)
%--------------------------------------------------
%  df-feb-2k
%	given set of vectors, calculate sampling tensor
%  vects = [x1,y1,z1;x2,y2,z2;...;xN,yN,zN];
%  kfull = 1 ==> full tensor; =0 ==> unique elements
%--------------------------------------------------
if nargin <2, kfull =1; end	%default = full tensor
if kfull ==1,
   SS=zeros(3);
   SS=diag((3*mean(vects.^2)-1)/2);
   SS(1,2)=mean(vects(:,1).*vects(:,2))*3/2;
   SS(1,3)=mean(vects(:,1).*vects(:,3))*3/2;
   SS(2,3)=mean(vects(:,2).*vects(:,3))*3/2;
   SS(3,2)=SS(2,3);
   SS(2,1)=SS(1,2);
   SS(3,1)=SS(1,3);
else
   SS=zeros(1,6);
   SS(1:3)=(3*mean(vects.^2)-1)/2;
   SS(4)=mean(vects(:,1).*vects(:,2))*3/2;
   SS(5)=mean(vects(:,1).*vects(:,3))*3/2;
   SS(6)=mean(vects(:,2).*vects(:,3))*3/2;
end
   
return 
%==================================================