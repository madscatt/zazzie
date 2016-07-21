function z=combine2(A,B,reslst,col,kNaN),
%----------------------------------
%  df-jan-98	df-dec-97
%	combine two arrays
%	according to the reslst
%	col -- data-column to combine
%	kNaN =1 to remove all NaN's
%----------------------------------
if nargin<5, kNaN=0;	end	%default: no NaN remove
if nargin<4, col=2; 	end	%default: col=2
if nargin<3, reslst=[]; end	%default: all resid.
if isempty(reslst), 
  rlst=[min([A(:,1);B(:,1)]):max([A(:,1);B(:,1)])]';
else
  rlst=reslst(:);
end
nres=length(rlst);
z=NaN*ones(nres,3);
z(:,1)=rlst;
for ii=1:nres,
  indA=find(A(:,1)==rlst(ii));
  if ~isempty(indA), z(ii,2)=A(indA,col); end
  indB=find(B(:,1)==rlst(ii));
  if ~isempty(indB), z(ii,3)=B(indB,col); end
end
if kNaN==1,	%all NaN's must out!
  ind=(~isnan(z(:,2))).*(~isnan(z(:,3)));
  z=z(find(ind),:);
end
return
%===================================
