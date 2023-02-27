function z=ftest(var1,n1,var2,n2)
%----------------------------------------------------------
%   df-nov-95   recovered
%
%       produces F-test: given variances, var1 and var2,
%       and numbers of degrees of freedom, n1 and n2,
%       calculates the probability that the
%       variances in the two data sets are actually
%       consistent  (the null-hypothesis):
%       if the probability is low, the null-hypothesis
%       is rejected (see Numer.Recipes, Ch.14.2)
%       INPUT:  var1,var2 - variances
%               n1,n2 - number of degr. of freedom
%
%			FOR F-TESTING one fit against the other, 
%			i.e. F=[(SSW-SSE)/(dfW-dfE)]/[SSE/dfE],
%				e.g. SSE is a better fit than SSW, SSE<SSW, dfW>dfE
%			input var1=(SSW-SSE)/(dfW-dfE)
%					n1=dfW-dfE (degr.of freedom of the numerator)
%					var2=SSE/dfE
%					n2=dfE	(degrees of freedom of the denominator)
%
% Copyright (c) 1996-97 by David Fushman, The Rockefeller University
%----------------------------------------------------------
if var1>var2,
    F=var1/var2;
    df1=n1;
    df2=n2;
else
    F=var2/var1;
    df1=n2;
    df2=n1;
end
%F
x=df2/(df2+df1*F);
z=betainc(x,df2/2,df1/2);       %incomplete BETA function
%it has to be multiplied by 2 for the two-tail distribution
return
%==========================================================                 