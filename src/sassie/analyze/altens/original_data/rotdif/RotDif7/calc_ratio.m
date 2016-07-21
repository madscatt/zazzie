function z=calc_ratio(rotmatrix,coord,wN,par)
%---------------------------------------------------------------------------
%
%----------------------------------------------------------------------------

Dxx=par(1)*1.0e07;
Dyy=par(2)*1.0e07;
Dzz=par(3)*1.0e07;



%------Define D--------------------------------------------------------------------

Diso=(Dxx+Dyy+Dzz)/3;
Dsq=(Dxx*Dyy+Dyy*Dzz+Dxx*Dzz)/3;

D1=4*Dxx+Dyy+Dzz;
D2=Dxx+4*Dyy+Dzz;
D3=Dxx+Dyy+4*Dzz;
D4=6*Diso+6*sqrt(Diso^2-Dsq);
D5=6*Diso-6*sqrt(Diso^2-Dsq);


coord_r=coord;
   
%-------di(i=x,y,z)----------------------------------------------------------------
   
dx=(Dxx-Diso)/sqrt(Diso^2-Dsq);
dy=(Dyy-Diso)/sqrt(Diso^2-Dsq);
dz=(Dzz-Diso)/sqrt(Diso^2-Dsq);

%-------Ai(i=1..5)-----------------------------------------------------------------

res1=(1/4)*(3*(coord_r(:,1).^4+coord_r(:,2).^4+coord_r(:,3).^4)-1);
res2=(1/12)*(dx*(3*coord_r(:,1).^4+6*(coord_r(:,2).^2).*(coord_r(:,3).^2)-1)...
      		+dy*(3*coord_r(:,2).^4+6*(coord_r(:,3).^2).*(coord_r(:,1).^2)-1)...
            +dz*(3*coord_r(:,3).^4+6*(coord_r(:,1).^2).*(coord_r(:,2).^2)-1));


A1=3*(coord_r(:,2).^2).*(coord_r(:,3).^2);
A2=3*(coord_r(:,1).^2).*(coord_r(:,3).^2); 
A3=3*(coord_r(:,1).^2).*(coord_r(:,2).^2); 
A4=res1-res2;
A5=res1+res2;
   

%-------Calculate the spectral density functions---------------------------------
   
j_zero=A1/D1+A2/D2+A3/D3+A4/D4+A5/D5;
j_wN=A1*(D1/(wN^2+D1^2))+A2*(D2/(wN^2+D2^2))+A3*(D3/(wN^2+D3^2))...
    +A4*(D4/(wN^2+D4^2))+A5*(D5/(wN^2+D5^2));
         
%-------Calculate the target function--------------------------------------------
   
z=(3/4)*(j_wN./j_zero);

   
return

%==================================================================================
