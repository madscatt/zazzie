function z=combine_all(reslst,NH,Di)
%---------------------------------------------------------------------
%	Combine residue list, res. dip. couplings and x,y,z co-ordinates 
%	
%
%	Format :
%	col 1 : residue number
%   col 2 : D res. dip. couplings
%	col 3 : x coordinate
%	col 4 : y coordinate
%	col 5 : z coordinates
%		
%		
%  The x,y,z co-ordinates are converted to unit vectors
% 	in the x,y and z directions.
%
% 	ow 07/02/2002
%---------------------------------------------------------------------


%reslst=reslst';
%reslist=reslst(:,1)
%pause
nres=length(reslst)
pause
z(:,1)=sort(reslst(:));


	for ii=2:5
      z(:,ii)=NaN*ones(1,nres)';
   end

   

for ii=1:nres   
   ii
   indcoor=find(NH(:,1)==reslst(ii));
   indDi=find(Di(:,1)==reslst(ii))
   pause
  
   
   if ~isempty(indcoor) & ~isempty(Di) 
      
      z(ii,2)=Di(indDi,2)
      pause
            
      
         
      	z(ii,3)=NH(indcoor,2)./sqrt(NH(indcoor,2).^2+NH(indcoor,3).^2+NH(indcoor,4).^2);
      	z(ii,4)=NH(indcoor,3)./sqrt(NH(indcoor,2).^2+NH(indcoor,3).^2+NH(indcoor,4).^2);
      	z(ii,5)=NH(indcoor,4)./sqrt(NH(indcoor,2).^2+NH(indcoor,3).^2+NH(indcoor,4).^2);
     
      
  	end 
     
 end

z(any(isnan(z)'),:)=[];

return

%=======================================================================
   
