function z=invcovar_iso2(tau_i,sigma,freq)


tau_i=tau_i*1e-5;

wN=freq*0.1013756*2*pi*1.0e06; 


dr_ci=((3/4)*(2*wN^2*tau_i))/(1+wN^2*tau_i^2)^2; 
z=dr_ci'*dr_ci;
   
return

%==================================================================================