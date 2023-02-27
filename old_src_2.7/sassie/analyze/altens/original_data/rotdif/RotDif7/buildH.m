function vHN=buildH(vCaN,vNC,angle),
%----------------------------------------------
%  df-may-98
%	given coordinates of the CaN and NC vectors
%	build the NH vector, which is in the peptide 
%	plane, at a ~122 deg angle to NCa
%	zz=Zax/norm_zz
%	xx=[REF x zz]/norm_xx
%	yy=[zz x xx]
%  NH = Rot of vv around xx by (109.8-90)
%  ALL input vectors are [1 x 3] - vectors!!!
%  vCaN is directed from N to Ca;
%  vNC is directed from C to N (!)
%  angle is angle between CaN and HN
%     C\  |xx
%	  \ |
% 	   \|--------->Ca
%	   /\N        zz
% 	  /  \
% HN /      \ 
%            _|yy REF
%---------------------------------------------
if nargin < 3, angle=121.9947;	end	%default INSIGHT
xx=zeros(1,3);				%space/size holder
yy=zeros(1,3);				%space/size holder
zz=zeros(1,3);				%space/size holder
nr1=size(vCaN,1);nr2=size(vNC,1);
vHN=zeros(nr1,3);
for ii=1:nr1,
   Zax=vCaN(ii,:);
   REF=vNC(ii,:);
   %--------------normalize the Z-axis------------
   norm=sqrt(Zax*Zax');
   zz=Zax/norm;
   %---------------create the X-axis--------------
   %xx=cross(REF,zz) should do the same job
   xx(1)=zz(3)*REF(2)-REF(3)*zz(2);  	% xx=[REF x zz]
   xx(2)=zz(1)*REF(3)-REF(1)*zz(3);
   xx(3)=zz(2)*REF(1)-REF(2)*zz(1);
   norm=sqrt(xx*xx');   		%normalize xx 
   xx=xx/norm;

   %---------------create the Y-axis--------------
   yy(1)=zz(2)*xx(3)-xx(2)*zz(3);     %yy = [zz x xx] 
   yy(2)=zz(3)*xx(1)-xx(3)*zz(1);
   yy(3)=zz(1)*xx(2)-xx(1)*zz(2);
   %rotation around xx
  % R=eye(3);
  %phi=-(122-90)/180*pi; %-(109.8-90)/180*pi;	      %rotation angle
  %R=[1 0 0;0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
  %vHN(ii,:)=(R*yy')';   
  rrr=rotate_df([xx;yy;zz],90-angle,1);	%121.9947
  vHN(ii,:)=rrr(2,:);
end   
return
%==============================================
%the following should be checked: xx*yy'=xx*zz'=zz*yy'=0
%                                 xx*xx'=yy*yy'=zz*zz'=1
