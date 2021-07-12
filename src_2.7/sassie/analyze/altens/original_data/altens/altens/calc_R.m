function R_coeff=calc_R(exp,calc)



Dobs_Dcalc=exp-calc;
num=Dobs_Dcalc.*Dobs_Dcalc;
num=mean(num);


denom=2*mean(exp.*exp);

R_coeff=sqrt(num/denom);

return


