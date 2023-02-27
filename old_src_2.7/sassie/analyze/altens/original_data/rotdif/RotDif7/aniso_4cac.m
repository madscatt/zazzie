function k=aniso(v1,v2,v3);

% calculates anisotropy from Dxx, Dyy, Dzz
%force to use prolate case:

    k=2*v3/(v1+v2);
return

if abs(v1-v2) <= abs(v2-v3),  % prolate case
    k=2*v3/(v1+v2);
else                          % oblate case
    k=2*v1/(v2+v3);
end
return