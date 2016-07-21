function z=matrix_dalpha(alpha,beta,gamma)

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

z=[-cg*cb*sa-sg*ca ca*cb*cg-sg*sa 0;sg*cb*sa-cg*ca -sg*cb*ca-cg*sa 0;-sb*sa sb*ca 0];

return

%=====================================================
