function mat=invcovar_fullax_angle(par,sigma,rotmatrix,dR_dalpha,dR_dbeta,dR_dgamma,coord,wN)

%	Returns the inverse covariance matrix for the fully
%	asymmetric case. Input parameters :
%	par(1) --> Dxx
%	par(2) --> Dyy
%	par(3) --> Dzz
%	rotmatrix --> rotation matrix at chi^2 minimum
%
%	ow-University of Maryland- 2002



Dxx=par(1)*1e7;
Dyy=par(2)*1e7;
Dzz=par(3)*1e7;



%------Define D--------------------------------------------------------------------

Diso=(Dxx+Dyy+Dzz)/3;
Dsq=(Dxx*Dyy+Dyy*Dzz+Dxx*Dzz)/3;
denom=(Diso^2-Dsq);

D1=4*Dxx+Dyy+Dzz;
D2=Dxx+4*Dyy+Dzz;
D3=Dxx+Dyy+4*Dzz;
D4=6*Diso+6*sqrt(Diso^2-Dsq);
D5=6*Diso-6*sqrt(Diso^2-Dsq);

%-------Rotate the co-ordinates to the diffusion frame-----------------------------

coord_r=(rotmatrix*(coord'))';
   
%-------di(i=x,y,z)----------------------------------------------------------------
   
dx=(Dxx-Diso)/sqrt(Diso^2-Dsq);
dy=(Dyy-Diso)/sqrt(Diso^2-Dsq);
dz=(Dzz-Diso)/sqrt(Diso^2-Dsq);

%-------Ai(i=1..5)-----------------------------------------------------------------

res1=(1/4)*(3*(coord_r(:,1).^4+coord_r(:,2).^4+coord_r(:,3).^4)-1);
res2a=(3*coord_r(:,1).^4+6*(coord_r(:,2).^2).*(coord_r(:,3).^2)-1);
res2b=(3*coord_r(:,2).^4+6*(coord_r(:,3).^2).*(coord_r(:,1).^2)-1);
res2c=(3*coord_r(:,3).^4+6*(coord_r(:,1).^2).*(coord_r(:,2).^2)-1);
res2=(1/12)*(dx*res2a+dy*res2b+dz*res2c);




A1=3*(coord_r(:,2).^2).*(coord_r(:,3).^2);
A2=3*(coord_r(:,1).^2).*(coord_r(:,3).^2); 
A3=3*(coord_r(:,1).^2).*(coord_r(:,2).^2); 
A4=res1-res2;
A5=res1+res2;
   

%-------Calculate the spectral density functions---------------------------------
   
j_zero=A1/D1+A2/D2+A3/D3+A4/D4+A5/D5;
j_wN=A1*(D1/(wN^2+D1^2))+A2*(D2/(wN^2+D2^2))+A3*(D3/(wN^2+D3^2))...
    +A4*(D4/(wN^2+D4^2))+A5*(D5/(wN^2+D5^2));
         



%==================================================================================
%            derivative of the spectral density components
%==================================================================================

d_coor_dalpha=(dR_dalpha*(coord'))';
d_coor_dbeta=(dR_dbeta*(coord'))';
d_coor_dgamma=(dR_dgamma*(coord'))';


dres1_dalpha=1/4*(3*(4*(coord_r(:,1).^3).*d_coor_dalpha(:,1)+4*(coord_r(:,2).^3).*d_coor_dalpha(:,2)+4*(coord_r(:,3).^3).*d_coor_dalpha(:,3))-1);
dres1_dbeta=1/4*(3*(4*(coord_r(:,1).^3).*d_coor_dbeta(:,1)+4*(coord_r(:,2).^3).*d_coor_dbeta(:,2)+4*(coord_r(:,3).^3).*d_coor_dbeta(:,3))-1);
dres1_dgamma=1/4*(3*(4*(coord_r(:,1).^3).*d_coor_dgamma(:,1)+4*(coord_r(:,2).^3).*d_coor_dgamma(:,2)+4*(coord_r(:,3).^3).*d_coor_dgamma(:,3))-1);

dres2a_dalpha=12*(coord_r(:,1).^3).*d_coor_dalpha(:,1)+6*(2*coord_r(:,2).*d_coor_dalpha(:,2).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dalpha(:,3).*(coord_r(:,2).^2));
dres2a_dbeta=12*(coord_r(:,1).^3).*d_coor_dbeta(:,1)+6*(2*coord_r(:,2).*d_coor_dbeta(:,2).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dbeta(:,3).*(coord_r(:,2).^2));
dres2a_dgamma=12*(coord_r(:,1).^3).*d_coor_dgamma(:,1)+6*(2*coord_r(:,2).*d_coor_dgamma(:,2).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dgamma(:,3).*(coord_r(:,2).^2));


dres2b_dalpha=12*(coord_r(:,2).^3).*d_coor_dalpha(:,2)+6*(2*coord_r(:,1).*d_coor_dalpha(:,1).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dalpha(:,3).*(coord_r(:,1).^2));
dres2b_dbeta=12*(coord_r(:,2).^3).*d_coor_dbeta(:,2)+6*(2*coord_r(:,1).*d_coor_dbeta(:,1).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dbeta(:,3).*(coord_r(:,1).^2));
dres2b_dgamma=12*(coord_r(:,2).^3).*d_coor_dgamma(:,2)+6*(2*coord_r(:,1).*d_coor_dgamma(:,1).*(coord_r(:,3).^2)+2*coord_r(:,3).*d_coor_dgamma(:,3).*(coord_r(:,1).^2));


dres2c_dalpha=12*(coord_r(:,3).^3).*d_coor_dalpha(:,3)+6*(2*coord_r(:,2).*d_coor_dalpha(:,2).*(coord_r(:,1).^2)+2*coord_r(:,1).*d_coor_dalpha(:,1).*(coord_r(:,2).^2));
dres2c_dbeta=12*(coord_r(:,3).^3).*d_coor_dbeta(:,3)+6*(2*coord_r(:,2).*d_coor_dbeta(:,2).*(coord_r(:,1).^2)+2*coord_r(:,1).*d_coor_dbeta(:,1).*(coord_r(:,2).^2));
dres2c_dgamma=12*(coord_r(:,3).^3).*d_coor_dgamma(:,3)+6*(2*coord_r(:,2).*d_coor_dgamma(:,2).*(coord_r(:,1).^2)+2*coord_r(:,1).*d_coor_dgamma(:,1).*(coord_r(:,2).^2));



dres2_dalpha=1/12*(dx*dres2a_dalpha+dy*dres2b_dalpha+dz*dres2c_dalpha);
dres2_dbeta=1/12*(dx*dres2a_dbeta+dy*dres2b_dbeta+dz*dres2c_dbeta);
dres2_dgamma=1/12*(dx*dres2a_dgamma+dy*dres2b_dgamma+dz*dres2c_dgamma);







dA1_dalpha=3*(2*coord_r(:,2).*(coord_r(:,3).^2).*d_coor_dalpha(:,2)+2*(coord_r(:,2).^2).*coord_r(:,3).*d_coor_dalpha(:,3));
dA1_dbeta=3*(2*coord_r(:,2).*(coord_r(:,3).^2).*d_coor_dbeta(:,2)+2*(coord_r(:,2).^2).*coord_r(:,3).*d_coor_dbeta(:,3));
dA1_dgamma=3*(2*coord_r(:,2).*(coord_r(:,3).^2).*d_coor_dgamma(:,2)+2*(coord_r(:,2).^2).*coord_r(:,3).*d_coor_dgamma(:,3));


dA2_dalpha=3*(2*coord_r(:,1).*(coord_r(:,3).^2).*d_coor_dalpha(:,1)+2*(coord_r(:,1).^2).*coord_r(:,3).*d_coor_dalpha(:,3));
dA2_dbeta=3*(2*coord_r(:,1).*(coord_r(:,3).^2).*d_coor_dbeta(:,1)+2*(coord_r(:,1).^2).*coord_r(:,3).*d_coor_dbeta(:,3));
dA2_dgamma=3*(2*coord_r(:,1).*(coord_r(:,3).^2).*d_coor_dgamma(:,1)+2*(coord_r(:,1).^2).*coord_r(:,3).*d_coor_dgamma(:,3));


dA3_dalpha=3*(2*coord_r(:,1).*(coord_r(:,2).^2).*d_coor_dalpha(:,1)+2*(coord_r(:,1).^2).*coord_r(:,2).*d_coor_dalpha(:,2));
dA3_dbeta=3*(2*coord_r(:,1).*(coord_r(:,2).^2).*d_coor_dbeta(:,1)+2*(coord_r(:,1).^2).*coord_r(:,2).*d_coor_dbeta(:,2));
dA3_dgamma=3*(2*coord_r(:,1).*(coord_r(:,2).^2).*d_coor_dgamma(:,1)+2*(coord_r(:,1).^2).*coord_r(:,2).*d_coor_dgamma(:,2));


dA4_dalpha=dres1_dalpha-dres2_dalpha;
dA4_dbeta=dres1_dbeta-dres2_dbeta;
dA4_dgamma=dres1_dgamma-dres2_dgamma;

dA5_dalpha=dres1_dalpha+dres2_dalpha;
dA5_dbeta=dres1_dbeta+dres2_dbeta;
dA5_dgamma=dres1_dgamma+dres2_dgamma;



dj0_dalpha=(1/(D1^2))*(D1*dA1_dalpha)+(1/(D2^2))*(D2*dA2_dalpha)+(1/(D3^2))*(D3*dA3_dalpha)+(1/(D4^2))*(D4*dA4_dalpha)+(1/(D5^2))*(D5*dA5_dalpha);
dj0_dbeta=(1/(D1^2))*(D1*dA1_dbeta)+(1/(D2^2))*(D2*dA2_dbeta)+(1/(D3^2))*(D3*dA3_dbeta)+(1/(D4^2))*(D4*dA4_dbeta)+(1/(D5^2))*(D5*dA5_dbeta);
dj0_dgamma=(1/(D1^2))*(D1*dA1_dgamma)+(1/(D2^2))*(D2*dA2_dgamma)+(1/(D3^2))*(D3*dA3_dgamma)+(1/(D4^2))*(D4*dA4_dgamma)+(1/(D5^2))*(D5*dA5_dgamma);


djw_dalpha=(D1/(D1^2+wN^2))*dA1_dalpha+(D2/(D2^2+wN^2))*dA2_dalpha+(D3/(D3^2+wN^2))*dA3_dalpha+(D4/(D4^2+wN^2))*dA4_dalpha+(D5/(D5^2+wN^2))*dA5_dalpha;
djw_dbeta=(D1/(D1^2+wN^2))*dA1_dbeta+(D2/(D2^2+wN^2))*dA2_dbeta+(D3/(D3^2+wN^2))*dA3_dbeta+(D4/(D4^2+wN^2))*dA4_dbeta+(D5/(D5^2+wN^2))*dA5_dbeta;
djw_dgamma=(D1/(D1^2+wN^2))*dA1_dgamma+(D2/(D2^2+wN^2))*dA2_dgamma+(D3/(D3^2+wN^2))*dA3_dgamma+(D4/(D4^2+wN^2))*dA4_dgamma+(D5/(D5^2+wN^2))*dA5_dgamma;



dr_dalpha=((3/4)*(j_zero.*(djw_dalpha)-j_wN.*(dj0_dalpha))./(j_zero).^2)./sigma;
dr_dbeta=((3/4)*(j_zero.*(djw_dbeta)-j_wN.*(dj0_dbeta))./(j_zero).^2)./sigma;
dr_dgamma=((3/4)*(j_zero.*(djw_dgamma)-j_wN.*(dj0_dgamma))./(j_zero).^2)./sigma;






dDiso_dDxx=1/3;
dDiso_dDyy=1/3;
dDiso_dDzz=1/3;


dDsq_dDxx=1/3*(Dyy+Dzz);
dDsq_dDyy=1/3*(Dxx+Dzz);
dDsq_dDzz=1/3*(Dyy+Dxx);


dx_dDxx=(2/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDxx)*(Dxx-Diso));
dx_dDyy=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDyy)*(Dxx-Diso));
dx_dDzz=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDzz)*(Dxx-Diso));

dy_dDxx=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDxx)*(Dyy-Diso));
dy_dDyy=(2/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDyy)*(Dyy-Diso));
dy_dDzz=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDzz)*(Dyy-Diso));

dz_dDxx=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDzz)*(Dzz-Diso));
dz_dDyy=(-1/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDyy)*(Dzz-Diso));
dz_dDzz=(2/3)/sqrt(denom)-(1/2)*(1/(Diso^2-Dsq)^(3/2))*((2/3*Diso-dDsq_dDzz)*(Dzz-Diso));



dD1_dDxx=4;
dD2_dDxx=1;
dD3_dDxx=1;
dD4_dDxx=2+(3/sqrt(denom))*(2/3*Diso-dDsq_dDxx);
dD5_dDxx=2-(3/sqrt(denom))*(2/3*Diso-dDsq_dDxx);



dD1_dDyy=1;
dD2_dDyy=4;
dD3_dDyy=1;
dD4_dDyy=2+(3/sqrt(denom))*(2/3*Diso-dDsq_dDyy);
dD5_dDyy=2-(3/sqrt(denom))*(2/3*Diso-dDsq_dDyy);

dD1_dDzz=1;
dD2_dDzz=1;
dD3_dDzz=4;
dD4_dDzz=2+(3/sqrt(denom))*(2/3*Diso-dDsq_dDzz);
dD5_dDzz=2-(3/sqrt(denom))*(2/3*Diso-dDsq_dDzz);


dA4_dDxx=(-1/12)*(res2a*dx_dDxx+res2b*dy_dDxx+res2c*dz_dDxx);
dA5_dDxx=(1/12)*(res2a*dx_dDxx+res2b*dy_dDxx+res2c*dz_dDxx);
dA4_dDyy=(-1/12)*(res2a*dx_dDyy+res2b*dy_dDyy+res2c*dz_dDyy);
dA5_dDyy=(1/12)*(res2a*dx_dDyy+res2b*dy_dDyy+res2c*dz_dDyy);
dA4_dDzz=(-1/12)*(res2a*dx_dDzz+res2b*dy_dDzz+res2c*dz_dDzz);
dA5_dDzz=(1/12)*(res2a*dx_dDzz+res2b*dy_dDzz+res2c*dz_dDzz);


dj0_dDxx=-(A1/D1^2)*dD1_dDxx-(A2/D2^2)*dD2_dDxx-(A3/D3^2)*dD3_dDxx+(1/D4^2)*(dA4_dDxx*D4-dD4_dDxx*A4)+(1/D5^2)*(dA5_dDxx*D5-dD5_dDxx*A5);
dj0_dDyy=-(A1/D1^2)*dD1_dDyy-(A2/D2^2)*dD2_dDyy-(A3/D3^2)*dD3_dDyy+(1/D4^2)*(dA4_dDyy*D4-dD4_dDyy*A4)+(1/D5^2)*(dA5_dDyy*D5-dD5_dDyy*A5);
dj0_dDzz=-(A1/D1^2)*dD1_dDzz-(A2/D2^2)*dD2_dDzz-(A3/D3^2)*dD3_dDzz+(1/D4^2)*(dA4_dDzz*D4-dD4_dDzz*A4)+(1/D5^2)*(dA5_dDzz*D5-dD5_dDzz*A5);


djw_dDxx=(A1/(D1^2+wN^2))*dD1_dDxx-dD1_dDxx*((2*D1^2*A1)/(D1^2+wN^2)^2)+(A2/(D2^2+wN^2))*dD2_dDxx-dD1_dDxx*((2*D2^2*A2)/(D2^2+wN^2)^2)...
    +(A3/(D3^2+wN^2))*dD3_dDxx-dD3_dDxx*((2*D3^2*A3)/(D3^2+wN^2)^2)+(D4/(D4^2+wN^2))*dA4_dDxx+(A4/(D4^2+wN^2))*dD4_dDxx...
    -((2*D4^2*A4)/(D4^2+wN^2)^2)*dD4_dDxx+(D5/(D5^2+wN^2))*dA5_dDxx+(A5/(D5^2+wN^2))*dD5_dDxx...
    -((2*D5^2*A5)/(D5^2+wN^2)^2)*dD5_dDxx;

djw_dDyy=(A1/(D1^2+wN^2))*dD1_dDyy-dD1_dDyy*((2*D1^2*A1)/(D1^2+wN^2)^2)+(A2/(D2^2+wN^2))*dD2_dDyy-dD1_dDyy*((2*D2^2*A2)/(D2^2+wN^2)^2)...
    +(A3/(D3^2+wN^2))*dD3_dDyy-dD3_dDyy*((2*D3^2*A3)/(D3^2+wN^2)^2)+(D4/(D4^2+wN^2))*dA4_dDyy+(A4/(D4^2+wN^2))*dD4_dDyy...
    -((2*D4^2*A4)/(D4^2+wN^2)^2)*dD4_dDyy+(D5/(D5^2+wN^2))*dA5_dDyy+(A5/(D5^2+wN^2))*dD5_dDyy...
    -((2*D5^2*A5)/(D5^2+wN^2)^2)*dD5_dDyy;

djw_dDzz=(A1/(D1^2+wN^2))*dD1_dDzz-dD1_dDzz*((2*D1^2*A1)/(D1^2+wN^2)^2)+(A2/(D2^2+wN^2))*dD2_dDzz-dD1_dDzz*((2*D2^2*A2)/(D2^2+wN^2)^2)...
    +(A3/(D3^2+wN^2))*dD3_dDzz-dD3_dDzz*((2*D3^2*A3)/(D3^2+wN^2)^2)+(D4/(D4^2+wN^2))*dA4_dDzz+(A4/(D4^2+wN^2))*dD4_dDzz...
    -((2*D4^2*A4)/(D4^2+wN^2)^2)*dD4_dDzz+(D5/(D5^2+wN^2))*dA5_dDzz+(A5/(D5^2+wN^2))*dD5_dDzz...
    -((2*D5^2*A5)/(D5^2+wN^2)^2)*dD5_dDzz;

%------- generate the elements of the Hessian by summing over the data points ------


dr_dDxx=((3/4)*(j_zero.*(djw_dDxx)-j_wN.*(dj0_dDxx))./(j_zero).^2)./sigma;
dr_dDyy=((3/4)*(j_zero.*(djw_dDyy)-j_wN.*(dj0_dDyy))./(j_zero).^2)./sigma;
dr_dDzz=((3/4)*(j_zero.*(djw_dDzz)-j_wN.*(dj0_dDzz))./(j_zero).^2)./sigma;


mat(1,1)=dr_dDxx'*dr_dDxx;
mat(1,2)=dr_dDxx'*dr_dDyy;
mat(1,3)=dr_dDxx'*dr_dDzz;
mat(1,4)=dr_dDxx'*dr_dalpha;
mat(1,5)=dr_dDxx'*dr_dbeta;
mat(1,6)=dr_dDxx'*dr_dgamma;

mat(2,1)=mat(1,2);
mat(2,2)=dr_dDyy'*dr_dDyy;
mat(2,3)=dr_dDyy'*dr_dDzz;
mat(2,4)=dr_dDyy'*dr_dalpha;
mat(2,5)=dr_dDyy'*dr_dbeta;
mat(2,6)=dr_dDyy'*dr_dgamma;


mat(3,1)=mat(1,3);
mat(3,2)=mat(2,3);
mat(3,3)=dr_dDzz'*dr_dDzz;
mat(3,4)=dr_dDzz'*dr_dalpha;
mat(3,5)=dr_dDzz'*dr_dbeta;
mat(3,6)=dr_dDzz'*dr_dgamma;

mat(4,1)=mat(1,4);
mat(4,2)=mat(2,4);
mat(4,3)=mat(3,4);
mat(4,4)=dr_dalpha'*dr_dalpha;
mat(4,5)=dr_dalpha'*dr_dbeta;
mat(4,6)=dr_dalpha'*dr_dgamma;


mat(5,1)=mat(1,5);
mat(5,2)=mat(2,5);
mat(5,3)=mat(3,5);
mat(5,4)=mat(4,5);
mat(5,5)=dr_dbeta'*dr_dbeta;
mat(5,6)=dr_dbeta'*dr_dgamma;


mat(6,1)=mat(1,6);
mat(6,2)=mat(2,6);
mat(6,3)=mat(3,6);
mat(6,4)=mat(4,6);
mat(6,5)=mat(5,6);
mat(6,6)=dr_dgamma'*dr_dgamma;



return

%==================================================================================
