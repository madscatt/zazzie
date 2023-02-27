function z=vect2cosmat(vect)
%-------------------------------------------------------------------
%	df-nov-2011
%	Convert vectors into matrix of direction cosines
%   vect = [X Y Z]
%
%---------------------------------------------------------------------

cosX = vect(:,1);
cosY = vect(:,2);
cosZ = vect(:,3);

z(:,1)=(cosY).^2-(cosX).^2; 					
z(:,2)=(cosZ).^2-(cosX).^2;					
z(:,3)=2*cosX.*cosY;				
z(:,4)=2*cosX.*cosZ;			
z(:,5)=2*cosY.*cosZ;

return

%=======================================================================
