function draw_vect(vcoor,vcolor,width)
%-----------------------------------------
%	df-oct-01
%	draw a given set of vectors as sticks 
%	in 3d, originating in origin
%	INPUT: vcoor= a nx3 matrix [x,y,z]
%			 vcolor= a 'str' variable
%-----------------------------------------
if nargin < 3, width = 1; end
if nargin < 2, vcolor = 'b';	end	%default
plot3([0 vcoor(1,1)],[0 vcoor(1,2)],[0 vcoor(1,3)],vcolor,'LineWidth',width); 
hold on
grid on
for ii=1:size(vcoor,1),
   plot3([0 vcoor(ii,1)],[0 vcoor(ii,2)],[0 vcoor(ii,3)],vcolor,'LineWidth',width);
end
%=========================================