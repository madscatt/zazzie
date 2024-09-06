function mat_D=gen_Di(A,B)

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


phi_X=A(:,1);
phi_Y=A(:,2);
phi_Z=A(:,3);



z(:,1)=(cos(phi_Y)).^2-(cos(phi_X)).^2; 					
z(:,2)=(cos(phi_Z)).^2-(cos(phi_X)).^2;					
z(:,3)=2*cos(phi_X).*cos(phi_Y);				
z(:,4)=2*cos(phi_X).*cos(phi_Z);			
z(:,5)=2*cos(phi_Y).*cos(phi_Z);





%---- generate Di matrix ----------
mat_D=z*B;





return

%=======================================================================
