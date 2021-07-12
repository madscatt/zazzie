function sel=select(atnam,at_res,reslst,atlst,atdim,offs),
%---------------------------------------------------
%  df-dec-97
%	select atoms according to atnam & reslst
%Copyright (c) 1998 by David Fushman, The Rockefeller University
%---------------------------------------------------
if nargin<6, offs=0; end		%default
atnampos=[7+offs:7+offs+atdim-1];
nat=size(atnam,1);
indsel=zeros(nat,1);
isel=1;
for ii=1:nat,
  if isempty(reslst),			%all resid
     if ~isempty(findstr(atnam(ii,atnampos),atlst)),
          indsel(isel)=ii; isel=isel+1;
     end
  else					%res. selection
   if ~isempty(find(at_res(ii,2)==reslst)),
     if isempty(atlst),			%all atoms
        indsel(isel)=ii; isel=isel+1;
     else				%sel.atoms
        if ~isempty(findstr(atnam(ii,atnampos),atlst)),
          indsel(isel)=ii; isel=isel+1;
        end
     end
   end
  end
end
sel=indsel(1:isel-1);
return
%===================================================

