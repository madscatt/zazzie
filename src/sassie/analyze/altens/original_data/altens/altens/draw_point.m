function draw_point(vcoor)
%-----------------------------------------
%	df-oct-01
%	draw a given set of vectors as sticks 
%	in 3d, originating in origin
%	INPUT: vcoor= a nx3 matrix [x,y,z]
%			 vcolor= a 'str' variable
%-----------------------------------------


hold on
grid on
for ii=1:size(vcoor,1),
   plot3(vcoor(ii,1),vcoor(ii,2),vcoor(ii,3),'r.')
end
%=========================================