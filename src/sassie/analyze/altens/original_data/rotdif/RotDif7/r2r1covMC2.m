function z=r2r1covMC2(par,ratio,sigma,vcoor,wN)
%---------------------------------------
%	ow-2002-University of Maryland
%   based on the program r2r1 of df
%	calculate 1/[covariance matrix]
%	using first derivatives of 
%	f(r2/r1) 
%---------------------------------------
kexpdat=~isnan(ratio(:));
s=sigma(:);
Tx=par(1);
ee=par(2)-1;
wT2=wN*wN*Tx*Tx;

phi=par(3);      %phi=param(3)
theta=par(4);    %theta=param(4)
rot=ones(1,3);
rot=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
costheta=(rot*vcoor')';
sin2=1-costheta.^2;

dcostheta_dtheta=[cos(theta)*cos(phi) sin(phi)*cos(theta) -sin(theta)];
dcostheta_dphi=[-sin(phi)*sin(theta) cos(phi)*sin(theta) 0];

drot_dtheta=(dcostheta_dtheta*vcoor')';
drot_dphi=(dcostheta_dphi*vcoor')';

A=(3/4)/(1+wT2);
B=wT2/(wT2+(1+(1/6)*ee)^2);
C=3+2*ee;
D=4+3*ee+(2/9)*ee^2;
E=1+((4+(11/3)*ee+(19/18)*ee*ee+(5/54)*ee*ee)/(wT2+(1+(2/3)*ee)^2));



dsin2_dtheta=-2*costheta.*drot_dtheta;
dsin2_dphi=-2*costheta.*drot_dphi;


num=ee*sin2.*(D-E*ee*sin2);
denom=C+(1+(1/3)*ee*(2-3*sin2)).^2;

dnum_dtheta=D*ee*dsin2_dtheta-2*E*ee*ee*sin2.*dsin2_dtheta;
dnum_dphi=D*ee*dsin2_dphi-2*E*ee*ee*sin2.*dsin2_dphi;
ddenom_dtheta=-2*ee*dsin2_dtheta-(4/3)*ee*ee*dsin2_dtheta+2*ee*sin2.*dsin2_dtheta;
ddenom_dphi=-2*ee*dsin2_dphi-(4/3)*ee*ee*dsin2_dphi+2*ee*sin2.*dsin2_dphi;




dYdtheta=((A*B)*(dnum_dtheta.*(denom)-ddenom_dtheta.*(num))./(denom).^2)./s;
dYdphi=((A*B)*(dnum_dphi.*(denom)-ddenom_dphi.*(num))./(denom).^2)./s;


dYdT=2*wN*wN*Tx*3/4/(1+wT2)*...
   (-1/(1+wT2)*...
   (1+ee*sin2.*wT2/(wT2+(1+ee/6)^2)./...
   (3+2*ee+(1+ee/3*(2-3*sin2)).^2).*...
   (4+3*ee+2/9*ee*ee-ee*sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/...
   (wT2+(1+2*ee/3)^2))))+...
   ee*sin2.*(1+ee/6)^2/((wT2+(1+ee/6)^2)^2)./...
   (3+2*ee+(1+ee/3*(2-3*sin2)).^2).*...
   (4+3*ee+2/9*ee*ee-ee*sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/...
   (wT2+(1+2*ee/3)^2)))+...
   ee^2*sin2.^2.*wT2/(wT2+(1+ee/6)^2)./(3+2*ee+(1+ee/3*(2-3*sin2)).^2).*...
   (1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/((wT2+(1+2*ee/3)^2)^2)))./s; 

der_e1=1/(wT2+(1+ee/6)^2)./(3+2*ee+(1+ee/3*(2-3*sin2)).^2)-...
   ee/3*(1+ee/6)/((wT2+(1+ee/6)^2)^2)./(3+2*ee+(1+ee/3*(2-3*sin2)).^2)-...
   ee/(wT2+(1+ee/6)^2)*2*(1+1/3*(2-3*sin2).*(1+ee/3*(2-3*sin2)))./...
   ((3+2*ee+(1+ee/3*(2-3*sin2)).^2).^2);

der_e2=3+4/9*ee-sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/(wT2+(1+2*ee/3)^2))-...
   ee*sin2*(11/3+19/9*ee+5/18*ee*ee)/(wT2+(1+2*ee/3)^2)-...
   ee*sin2*(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)*4/3*(1+2*ee/3)/((wT2+(1+2*ee/3)^2)^2);

dYde=3/4/(1+wT2)*sin2.*wT2.*(der_e1.*...
  (4+3*ee+2/9*ee*ee-ee*sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/(wT2+(1+2*ee/3)^2)))+...
  ee/(wT2+(1+ee/6)^2)./(3+2*ee+(1+ee/3*(2-3*sin2)).^2).*der_e2)./s;


J0=[dYdT,dYde,dYdtheta,dYdphi];
J=J0(find(kexpdat),:);
z=J'*J; 			%inv.covar.matrix
return
%=================================================  












