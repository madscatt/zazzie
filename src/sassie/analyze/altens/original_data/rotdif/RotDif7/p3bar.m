function p3bar(x,y,z,h,siz,clr)
% function p3bar(x,y,z,h,siz,clr)
% P3BAR puts a 3d rectangular bar on your plot at position (x,y,z)
% height h  with a base square centered on x,y,z of size siz x siz
% color is clr in either [r,g,b] or in the one letter code used in plot
%  (for example: 'r' = red or [1 0 0] = red)
%
% There may be a better way to create the bar objects, but this worked
% for me.
%
% added enhancement to use vector inputs for x,y,z,&h and optionally clr
% 4/26/94  DKB (Dan Braithwaite)
%
%  EXAMPLE: x = [1:3]';
%y = [1:3]';
%z = zeros(3,1);
%h = [1:3]';
%clr = ['r';'g';'m'];
%p3bar(x,y,z,h,.1,clr)
%
%The above example should produce a simple 3 column graph
%
%  NOTE:  it is a good idea to add a mesh to show a base z value (say 0)
%
%   written by DKB (dank)  11/2/92 U. of A. Dept. of HWR
if size(x,1) < size(x,2), x = x'; end
if size(y,1) < size(y,2), y = y'; end
if size(z,1) < size(z,2), z = z'; end
if size(h,1) < size(h,2), h = h'; end
if ischar(clr)
if size(clr,1) < size(clr,2)
clr = clr';
elseif size(clr,1) == 1
clr = clr(ones(1,size(x,1)),:);
end
elseif size(clr,1) == 1 & size(clr,2) == 3
clr = clr(ones(1,size(x,1)),:);
end 
hold on
for j=1:size(x,1)
  xt = [x(j)-(siz/2) x(j)+(siz/2) x(j)+(siz/2) x(j)-(siz/2)];
  yt = [y(j)-(siz/2) y(j)-(siz/2) y(j)+(siz/2) y(j)+(siz/2)];
  xb = [x(j)-(siz/2) x(j)+(siz/2) x(j)+(siz/2) x(j)-(siz/2)];
  yb = [y(j)-(siz/2) y(j)-(siz/2) y(j)+(siz/2) y(j)+(siz/2)];
  zb = [z(j) z(j) z(j) z(j)];
  zt = zb+h(j);
  zs = [z(j) z(j) z(j)+h(j) z(j)+h(j)];
  y1 = [y(j)-(siz/2),y(j)-(siz/2),y(j)-(siz/2),y(j)-(siz/2)];
  y2 = [y(j)+(siz/2),y(j)+(siz/2),y(j)+(siz/2),y(j)+(siz/2)];
  x1 = [x(j)-(siz/2),x(j)-(siz/2),x(j)-(siz/2),x(j)-(siz/2)];
  x2 = [x(j)+(siz/2),x(j)+(siz/2),x(j)+(siz/2),x(j)+(siz/2)];
  ys = [y(j)-(siz/2) y(j)+(siz/2) y(j)+(siz/2) y(j)-(siz/2)];
 
  fill3(xt,yt,zt,clr(j));
  fill3(xb,yb,zb,clr(j));
 
  fill3(xt,y1,zs,clr(j));
  fill3(xt,y2,zs,clr(j));
  fill3(x1,ys,zs,clr(j));
  fill3(x2,ys,zs,clr(j));
end