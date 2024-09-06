function h=rhomb(v1,v2,v3)

if abs(v1-v2) <= abs(v2-v3),  % prolate case
   num=v2-v1;
   den=(v3-(v1+v2)/2);
else                          % oblate case
   num=v2-v3;
   den=(v1-(v2+v3)/2);
end
h=3/2*num/den;
return