function [P,F]=ftestcmp(var1,df1,var2,df2)
%----------------------------------------------------------
%   df-oct-98   modif df-mar-04
%
%       produces F-test-based comparison of two models
%		  with different numbers of degrees of freedom, to
%		  assess significance of the improvement of fit
%		  given variances, var1 and var2, and numbers of df,
%		  df1, and df2, 
%		  calculates F=[(SSW-SSE)/(dfW-dfE)]/[SSE/dfE]
%		  e.g. SSW=var1,SSE=var2, dfW=n1, dfE=n2, if (var1>var2, n1>n2)
%       (see Draper & Smith, p. 105; dfE=d.f. for denominator
%		    						dfW-dfE =d.f. for numerator)
%		  and returns probability that this value of F
%		  could have occurred by chance
%       INPUT:  var1,var2 - variances
%               df1,df2 - number of degr. of freedom
%
% Copyright (c) 1996-97 by David Fushman, The Rockefeller University
%----------------------------------------------------------
P=NaN; F=NaN;
if var1*df1>var2*df2,
   SSW=var1*df1; SSE=var2*df2; dfW=df1; dfE=df2; 
elseif var1*df1<var2*df2,
   SSW=var2*df2; SSE=var1*df1; dfW=df2; dfE=df1;
else								%var1*df1=var2*df2
   disp('equal target functions. cant perform the test!!!');
   return
end
dfNUM=dfW-dfE;
if dfNUM == 0, disp('zero d.f.NUM. can"t perform the F-test!!!'); return; end
F=(SSW-SSE)/(dfW-dfE)/(SSE/dfE);
if dfNUM < 0, disp('F-test: model with fewer degrees of freedom has higher chi^2 --> to be rejected!'); return; end
P=ftest((SSW-SSE)/dfNUM,dfNUM,SSE/dfE,dfE);
if F < 1, P = 1 - P; end        %DF:Mar04
return
%==========