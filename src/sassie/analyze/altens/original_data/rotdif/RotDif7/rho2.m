function z=rho2(param,ratio,sigma,vcoor,wN)

% function involved in the fitting procedure of
% axially symmetric case for LVM algorithm
% ow-2002-University of Maryland

Tx=param(1);                %tau X=param(1)
ee=param(2)-1;				%Dz2Dx=param(2);ee=Dz2Dx-1;
wT2=wN*wN*Tx*Tx;
phi=param(3);               %phi=param(3)
theta=param(4);             %theta=param(4)
rot=ones(1,3);
rot=[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
costheta=(rot*vcoor')';
sin2=1-costheta.^2;

%------------calculate ratio ---------
calc=3/4/(1+wT2)*(1+ee*sin2.*wT2/(wT2+(1+ee/6)^2)./...
 (3+2*ee+(1+ee/3*(2-3*sin2)).^2).*...
  (4+3*ee+2/9*ee*ee-ee*sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/...
  (wT2+(1+2*ee/3)^2))));



%------------get difference-----------
diff=(ratio-calc)./sigma;
z=diff;				%might need looping for MCC
return
%===============================================