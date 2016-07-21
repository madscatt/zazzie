function ff=ss2ff(ss)
%-------------------------------------
%  df-feb-2k
%   convert sampling vector/tensor into
%   fraction numbers
%-------------------------------------
[nr,nc]=size(ss);
if nr==1 & nc==6,		%vector!!!
	[stens,sdiag]=vec2tens(ss);
elseif nr==3 & nc==3, %tensor
   sdiag=eig(ss);
else
   error('wrong input!');
end
ff=(sdiag*2+1)/3;
return
%=====================================