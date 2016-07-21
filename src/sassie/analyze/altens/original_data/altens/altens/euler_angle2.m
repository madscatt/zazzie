function z=euler_angle2(r)


z(3)=atan2(r(2,3),r(1,3));
z(2)=atan2(sqrt(r(1,3)^2 + r(2,3)^2),r(3,3));
z(1)=atan2(r(3,2),-r(3,1));



z(1)=z(1)*180/pi;
z(2)=z(2)*180/pi;
z(3)=z(3)*180/pi;



return
