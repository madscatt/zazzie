function z=invcovar_iso(d_i,sigma,freq)

d_iso=d_i*1e7;
%wN=freq*(2710.5/26750)*2*pi/1000;
wN=freq*0.1013756*2*pi*1.0e06;

   j_zero=1/(6*d_iso);
   j_wN=(6*d_iso)/(36*d_iso^2+wN^2);
   
   dj0_ci=-1/(6*d_iso^2);
   djw_ci=(6*(wN^2-36*d_iso^2))/((36*d_iso^2+wN^2)^2);
   
   dr_ci=((3/4)*(j_zero.*(djw_ci)-j_wN.*(dj0_ci))./(j_zero).^2)./sigma;
   
   z=dr_ci'*dr_ci;
   
return

%==================================================================================