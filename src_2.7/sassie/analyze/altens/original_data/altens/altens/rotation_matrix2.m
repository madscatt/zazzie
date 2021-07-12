function z=rotation_matrix(alpha,beta,gamma)

%----------------------------------------------------
%
%  Calculates the rotation matrix for a given set
%  of Euler angles alpha, beta and gamma. [IN DEGREES]
%   
%----------------------------------------------------

ar = alpha*pi/180;
br = beta*pi/180;
gr = gamma*pi/180;

%------Build matrix----------------------------------

z(1,1)=cos(ar)*cos(br)*cos(gr)-sin(ar)*sin(gr);
z(1,2)=sin(ar)*cos(br)*cos(gr)+cos(ar)*sin(gr);
z(1,3)=-sin(br)*cos(gr);
z(2,1)=-cos(ar)*cos(br)*sin(gr)-sin(ar)*cos(gr);
z(2,2)=-sin(ar)*cos(br)*sin(gr)+cos(ar)*cos(gr);
z(2,3)=sin(br)*sin(gr);
z(3,1)=cos(ar)*sin(br);
z(3,2)=sin(ar)*sin(br);
z(3,3)=cos(br);

return

%=====================================================
