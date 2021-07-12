function z=rotation_matrix(alpha,beta,gamma)

%----------------------------------------------------
%
%  Calculates the rotation matrix for a given set
%  of Euler angles alpha, beta and gamma.
%   
%----------------------------------------------------

ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);
cg=cos(gamma);
sg=sin(gamma);

%------Build matrix----------------------------------

La=[ca sa 0;-sa ca 0;0 0 1];
Lb=[cb 0 -sb;0 1 0;sb 0 cb];
Lg=[cg sg 0;-sg cg 0;0 0 1];
z=Lg*Lb*La;

return

%=====================================================
