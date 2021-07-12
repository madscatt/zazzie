function z=matrix_dbeta(alpha,beta,gamma)

%----------------------------------------------------
%
%  Derive the rotation matrix with respect to alpha
%  ow-University Of Maryland-2002
%   
%----------------------------------------------------

ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);
cg=cos(gamma);
sg=sin(gamma);

%------Build matrix----------------------------------

z=[-cg*sb*ca -cg*sb*sa -cg*cb;sg*sb*ca sg*sb*sa sg*cb;cb*ca cb*sa -sb];

return

%=====================================================
