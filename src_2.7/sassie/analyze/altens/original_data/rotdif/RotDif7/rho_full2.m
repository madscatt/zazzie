function f=rho_full2(x,ratio,sigma,vcoor,wN)

% function used to compute difference between exp and th data
% used with the LVM algorithm

alpha=x(4);
beta=x(5);
gamma=x(6);

rotmatrix=rotation_matrix(alpha,beta,gamma);
f=target_full4(ratio,sigma,rotmatrix,vcoor,wN,x);

return