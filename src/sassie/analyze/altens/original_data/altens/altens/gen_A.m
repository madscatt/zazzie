function [z,z_inv,w,v,con_num]=gen_A(angles)

%-------------------------------------------------------------------
%	
%	Calculates the A matrix from the input co-ordinate
%	set and obtains the singular value decomposition
%	A = U.W.transpose(V). The elements of W are the 
%	singular values of A. The inverse of A is obtained
%	as inverse(A)=V.diagonal(1/W).transpose(U).
%
%	ow 07/02/2002
%
%---------------------------------------------------------------------

TOL=0.01;
phi_X=angles(:,1);
phi_Y=angles(:,2);
phi_Z=angles(:,3);




z(:,1)=(cos(phi_Y)).^2-(cos(phi_X)).^2; 					
z(:,2)=(cos(phi_Z)).^2-(cos(phi_X)).^2;					
z(:,3)=2*cos(phi_X).*cos(phi_Y);				
z(:,4)=2*cos(phi_X).*cos(phi_Z);			
z(:,5)=2*cos(phi_Y).*cos(phi_Z);


con_num=cond(z);					% Condition number

[u,w,v]=svd(z,0);					% z=u.w.transpose(v)

max_w=max(max(w));


for i=1:5
   if w(i,i) > TOL*max_w;
      s(i,i)=1/(w(i,i));				% Inverse of singular values
else 
      s(i,i)=0.0;					% Exclude very small values								
end      
end
z_inv=v*s*transpose(u);					% v.s.transpose(u)

return

%=======================================================================
