function orts1=rotate(orts,angle,axis)
%----------------------------------------------
%  df-feb-2k
%  rotation by an angle <angle> around one
%  of the axes <axis> of an arbitrary 
%  coord.system <orts>
%----------------------------------------------
cs=cos(angle*pi/180);
si=sin(angle*pi/180);
%determine rotation from <orts> to <xyx> frame:
R=orts;
%define rotation matrix aroun a selected axis:
if axis ==1, 		
    Rot=[ 1 0 0; 0 cs si; 0 si -cs]; 	%(around x) CHECK!!!   
elseif axis == 2, 	
    Rot=[ cs 0 -si; 0 1 0; si 0 cs]; 	%(around y)   
else			
    Rot=[ cs si 0;-si cs 0; 0 0 1];    %(around z) 
end
orts1=orts*R'*Rot*R;
%orts1=(R'*Rot*R*orts')'; this works as well
return
%==============================================

