function z=tau_ani_rho2Ds(tc,ani,rho)
%------------------------------------------------
%   df-apr-13;  df-aug-03
%   given TAUc, anisotropy and rhombicity,
%   calculate the principal values of the diff.tensor
%------------------------------------------------
if ani >= 1,                  %prolate model
    Dz=ani/2/TAU/(2+ani);
    Dx=(1-rho*(ani-1)/3)/2/TAU/(2+ani);
    Dy=(1+rho*(ani-1)/3)/2/TAU/(2+ani);
else                          %oblate model
    Dx=ani/2/TAU/(2+ani);
    Dz=(1-rho*(ani-1)/3)/2/TAU/(2+ani);
    Dy=(1+rho*(ani-1)/3)/2/TAU/(2+ani);
end
z=[Dx,Dy,Dz]*1e9
%================================================