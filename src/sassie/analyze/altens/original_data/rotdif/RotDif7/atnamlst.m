function [atlst,atdim]=atnamlst(atnam)
%-----------------------------------------
%  df-dec-97
%	make  a list of atom names
%
%Copyright (c) 1997 by David Fushman, The Rockefeller University
%-----------------------------------------
if ~isempty(atnam),                     %at-list options
   if strcmp(atnam,'bb'),               %backbone heavy
        atlst=['N CAC O S '];atdim=2;
   elseif strcmp(atnam,'hv')|strcmp(atnam,'heavy'),
        atlst=['NCOS']; atdim=1;
   else                                 %convert to line
       [natlst,atdim]=size(atnam);
       atlst=setstr(32*ones(1,natlst*atdim));
       for ii=1:natlst,
         atlst((ii-1)*atdim+1:ii*atdim)=atnam(ii,:);
       end
   end
else
  atdim=0;
end
return
%==========================================
