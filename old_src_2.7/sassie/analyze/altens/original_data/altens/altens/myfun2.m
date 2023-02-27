function F=myfun2(x,C);


ar=x(1,1);
br=x(1,2);
gr=x(1,3);




calc(1,1)=cos(ar)*cos(br)*cos(gr)-sin(ar)*sin(gr);
calc(1,2)=sin(ar)*cos(br)*cos(gr)+cos(ar)*sin(gr);
calc(1,3)=-sin(br)*cos(gr);
calc(1,4)=-cos(ar)*cos(br)*sin(gr)-sin(ar)*cos(gr);
calc(1,5)=-sin(ar)*cos(br)*sin(gr)+cos(ar)*cos(gr);
calc(1,6)=sin(br)*sin(gr);
calc(1,7)=cos(ar)*sin(br);
calc(1,8)=sin(ar)*sin(br);
calc(1,9)=cos(br);


z(1,1)=C(1,1);
z(1,2)=C(1,2);
z(1,3)=C(1,3);
z(1,4)=C(2,1);
z(1,5)=C(2,2);
z(1,6)=C(2,3);
z(1,7)=C(3,1);
z(1,8)=C(3,2);
z(1,9)=C(3,3);

diff=(z-calc);

F=diff;
return
