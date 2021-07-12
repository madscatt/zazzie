function[V] = atamult(A,Y,flag,varargin);
%ATAMULT Example Jacobian-matrix multiply
%
%	V = ATAMULT(Y,A) computes V = (A'*(A*Y)).
%
%	V = ATAMULT(Y,A,dummy,flag) computes V = (A'*(A*Y)) if flag = 0,
%                                            V = A*Y if flag > 0,
%                                            V =  A'*Y if flag < 0.
%
% Note: varargin is not used but must be provided in case 
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

%   Copyright 1990-2000 The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 2000/06/16 22:12:15 $

if nargin < 3 | flag == 0
   V = (A'*(A*Y));
elseif flag > 0
   V = A*Y;
else
   V = A'*Y;
end


