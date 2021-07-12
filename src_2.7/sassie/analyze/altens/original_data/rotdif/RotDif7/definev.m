function [v,dv]= definev(g,x,l,u);
%DEFINEV Scaling vector and derivative
%
%	[v,dv]= DEFINEV(g,x,l,u) returns v, distances to the
%   bounds corresponding to the sign of the gradient g, where
%   l is the vector of lower bounds, u is the vector of upper 
%   bounds. Vector dv is 0-1 sign vector (See ?? for more detail.)
%

%   Copyright 1990-2000 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2000/06/16 22:12:21 $

n = length(x); 
v = zeros(n,1); 
dv=zeros(n,1);
arg1 = (g < 0)  & (u <  inf ); 
arg2 = (g >= 0) & (l > -inf);
arg3 = (g < 0)  & (u == inf); 
arg4 = (g >= 0) & (l == -inf);
v(arg1)  = (x(arg1) - u(arg1)); 
dv(arg1) = 1;
v(arg2)  = (x(arg2) - l(arg2)); 
dv(arg2) = 1;
v(arg3)  = -1;
dv(arg3) = 0;
v(arg4)  = 1;
dv(arg4) = 0;

