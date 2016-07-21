function v=eulangMODIFIED(R)
% assumes R is in format [Rx;Ry;Rz] where Ri = [Rix,Riy,Riz]
% gets euler angles

Rt = R';
Dx=Rt(1:3,1);
Dy=Rt(1:3,2);
Dz=Rt(1:3,3);

beta = acos(Dz(3));
n = whichoctant(Dz);

if beta==0
    alpha = acos(Dx(1));
else
    alpha = acos(Dz(1)/sin(beta));
    if (1-abs(Dz(1)/sin(beta)))<10e-005
        alpha = real(acos(Dz(1)/sin(beta)));end
    if Dz(1)/sin(beta) > 0 & (n==6 | n==8)
        alpha = 2*pi - alpha;end
    if Dz(1)/sin(beta) < 0 & (n==5 | n==7)
        alpha = alpha + 2*(pi-alpha);end
end

xdubprime = [cos(alpha)*cos(beta);sin(alpha)*cos(beta);-sin(beta)];
a = dot(xdubprime,Dx);
b = dot(xdubprime,Dy);
gamma = acos(a);
if a>0 & b>0
    gamma = 2*pi - gamma;end
if a<0 & b>0
    gamma = gamma+2*(pi-gamma);end

if beta==0
    alpha = alpha+gamma;
    gamma = 0;
end
alpha = alpha*180/pi;
beta = beta*180/pi;
gamma = gamma*180/pi;
v = [alpha beta gamma];

return