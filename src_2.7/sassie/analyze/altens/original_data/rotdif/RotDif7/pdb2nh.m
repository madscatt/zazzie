function vNH=pdb2nh(fname,reslst,model)
%------------------------------------------------
%	df-oct-01
%	create NH-vectors 
%------------------------------------------------
if nargin < 3, model = 1; end		%default: 1st structure
if nargin < 2, reslst = []; end	%defaul: all resid

[coor,atnam,at_res]=readpdb(fname,reslst,[],model);
%--------first try the 'usual' naming convention
vNH=getNHvect(coor,atnam,at_res,reslst,0);
if ~isempty(vNH),
   return
end
   %---------try Insight
vNH=getNHvect(coor,atnam,at_res,reslst,-1);
if ~isempty(vNH),
   return
end
%---------build it yourself
disp('no amide Hs found in PDB, building NH vectors');
   
iN=select_at(atnam,at_res,reslst,' N ',3,-1);
iC=select_at(atnam,at_res,reslst,' C ',3,-1);
iCa=select_at(atnam,at_res,reslst,' CA',3,-1);
nN=size(iN,1);
vNH=NaN*ones(nN,4);
vNH(:,1)=at_res(iN,2);
angleCaNH=119.1; %121.9947;
vNHlen=1;
for ii=1:nN,
  if ~strcmp(atnam(iN(ii),11:13),'PRO'), %Non-PRO only!!!
     indC=find(at_res(iC,2)==at_res(iN(ii),2)-1);
     if ~isempty(indC),	%check if not N-term!		
         vCaN=coor(iCa(ii),:)-coor(iN(ii),:);
         vNC=coor(iN(ii),:)-coor(iC(indC),:);
         vNH(ii,2:4)=buildH(vCaN,vNC,angleCaNH)*vNHlen;
     end
  end
end
%--- weed NaNs-----
inotNaN=find(~isnan(vNH(:,2)));
vNH=vNH(inotNaN,:);
return
%===========================================