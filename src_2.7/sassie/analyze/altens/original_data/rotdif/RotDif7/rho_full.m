function f=rho_full(x,ratio,sigma,vcoor,wN)

%function used to calculate error in the fully anisotropic case
%ow-2002-University of Maryland

alpha=x(4);
beta=x(5);
gamma=x(6);

rotmatrix=rotation_matrix(alpha,beta,gamma);
f=target_full3(ratio,sigma,rotmatrix,vcoor,wN,x);

return