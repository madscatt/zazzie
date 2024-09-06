function d=whichoctant(D)
%-returns 1-8 for 3D vector organized [Dx;Dy;Dz]

x = D(1);
y = D(2);
z = D(3);

if x==0 | y==0 | z==0
    n=0;end

if x>0 & y>0 & z>0
    n=1;end
if x<0 & y>0 & z>0
    n=2;end
if x>0 & y>0 & z<0
    n=3;end
if x<0 & y>0 & z<0
    n=4;end
if x<0 & y<0 & z>0
    n=5;end
if x>0 & y<0 & z>0
    n=6;end
if x<0 & y<0 & z<0
    n=7;end
if x>0 & y<0 & z<0
    n=8;end

d=n;