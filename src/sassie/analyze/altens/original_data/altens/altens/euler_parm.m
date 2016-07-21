function [angle,zflag]=euler_parm(rin)

%===================================================================
%
%  Finds the Euler angles corresponding to the input rotation
%  matrix. REMEMBER that any row multiplied by -1 is also a
%  valid solution.
%
%  RG 04/23/2000
%===================================================================

TOL = 1.0e-12;
result=[];
angle=zeros(1,3);
zflag=0;
negative_flag=1;

if negative_flag > 0
	if det(rin) < 0
   	r=[rin(1,:);rin(2,:);-1*rin(3,:)];
	else
   	r=rin;
	end
else
   r=rin;
end

t=acos(r(3,3));

if sin(t) ~=0
   f1=asin(r(3,2)/sin(t));
   f2=pi-f1;
   k1=asin(r(2,3)/sin(t));
   k2=pi-k1;
   
   flag=0;
   
   flag=choose_euler(r,flag,f1,t,k1,1);
   flag=choose_euler(r,flag,f1,t,k2,2);
   flag=choose_euler(r,flag,f2,t,k1,3);
   flag=choose_euler(r,flag,f2,t,k2,4);
   
   result=[1,f1,t,k1;2,f1,t,k2;3,f2,t,k1;4,f2,t,k2];
   index=find(flag==result(:,1));
   if ~isempty(index)
      angle(1)=180*result(index,2)/pi;
      angle(2)=180*result(index,3)/pi;
      angle(3)=180*result(index,4)/pi;
   end
   for i=1:3
      if angle(i) < 0
         angle(i)=360-abs(angle(i));
      end
   end
   else
   if cos(t)+1 < TOL
      angle(2)=180;
     	angle(1)=(acos(-r(1,1)))*180/pi;
      angle(3)=0;
   end     
   if cos(t)-1 < TOL
      angle(2)=0;
      angle(1)=(acos(r(1,1)))*180/pi;   
      angle(3)=0;
   end
end

if angle(1)~=0 & angle(2)~=0 & angle(3)~=0  
   zflag=1;
end

return

%==========================================================================
