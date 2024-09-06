function z=ff2ksi(ff)
%--------------------------------------------
% df-mar-2k
%   given ff-set, calculate the generalized
%   sampling parameter KSI
%--------------------------------------------
[nr,nl]=size(ff);
if nl==4, ff=ff(:,2:4); end
z=((ff(:,1).^2+ff(:,2).^2+ff(:,3).^2)*3-1)/2;
return

   