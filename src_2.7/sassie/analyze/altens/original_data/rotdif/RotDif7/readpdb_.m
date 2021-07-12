function [coor,atnam,at_res]=readpdb(fname,reslst,atlst,model)
%-----------------------------------------------------------------
%   df-nov-97
%	read in a pdb data set 
%       uses: readasci.m
%
%Copyright (c) 1997 by David Fushman, The Rockefeller University
%-----------------------------------------------------------------
if nargin<4, model=1; 	end		%default: first model
if nargin<3, atlst=[];  end             %default: all atoms
if nargin<2, reslst=[]; end		%default: all resid.
pdb=readasci(fname);
disp(['got the ',fname,' data set, analyzing...']);
nlin=size(pdb,1);
select_at=zeros(nlin,1);
isel=1;
nmod=1;
termflag=0;				%term flag off
for ii=1:nlin,
  if strcmp('TER',pdb(ii,1:3))|strcmp('END',pdb(ii,1:3))|strcmp('.',pdb(ii,1)), 
     if termflag==0,
        nmod=nmod+1;
        termflag=1;			%term flag on
        if nmod>model, break; end
     end
  else
     termflag=0;			%term flag off
     if nmod==model			%read in the structure
       if strcmp('ATOM',pdb(ii,1:4)), 
         select_at(isel)=ii; isel=isel+1;
       end
     end
  end
end
sel=select_at(1:isel-1);
nat=length(sel);
if nat==0, error('no atoms found!!! wrong filename or model'); end  
disp([num2str(nat),' atoms read in']);
atnam=pdb(sel,8:26);
coor=zeros(nat,3);
at_res=zeros(nat,2);
at_res(:,1)=str2num(pdb(sel,7:11));
at_res(:,2)=str2num(pdb(sel,23:26));
coor(:,1:3)=str2num(pdb(sel,32:54));
return
%=================================================================
