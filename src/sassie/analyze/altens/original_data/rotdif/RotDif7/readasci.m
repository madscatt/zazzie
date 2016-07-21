function z=readasci(filename)
%---------------------------------------------------
%  df-feb-08 (took care of EOF problems) df-nov-97
%    read in the an ASCII file using fread
%---------------------------------------------------
fid = fopen(filename,'r');
F = fread(fid);
fclose(fid);
ENDLINE=10;
endl=(find(F==ENDLINE));
tot_len=length(F);
if endl(end)~=tot_len, endl=[endl;tot_len]; end  %if no endline at EOF
nlin=length(endl);
maxlen=endl(1);
for ii=2:nlin,
  if maxlen<endl(ii)-endl(ii-1), maxlen=endl(ii)-endl(ii-1);end
end
mtrx=ones(nlin,maxlen)*32;
mtrx(1,1:endl(1)-1)=F(1:endl(1)-1)';	%first line
for ii=2:nlin,
  mtrx(ii,1:endl(ii)-endl(ii-1)-1)=F(endl(ii-1)+1:endl(ii)-1)';
end
z=setstr(mtrx);
return
%====================================================
