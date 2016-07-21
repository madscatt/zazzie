function z=matrix_dgamma(alpha,beta,gamma)

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

z=[-sg*cb*ca-cg*sa -sg*cb*sa+cg*ca sg*sb;-cg*cb*ca+sg*sa -cg*cb*sa-sg*ca cg*sb;0 0 0];

return

%=====================================================
