function z=extract_par(S)
%--------------------------------------------
%   df-nov-2011
%   convert alignment tensor S into vector
%--------------------------------------------
z=NaN*zeros(5,1);
z(1) = S(2,2);
z(2) = S(3,3);
z(3) = S(1,2);
z(4) = S(1,3);
z(5) = S(2,3);

return