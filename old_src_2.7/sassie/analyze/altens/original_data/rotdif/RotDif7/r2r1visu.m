function [z,fitline,thetacrv,l_lim,u_lim,Q,corr,ratio_calc]=r2r1visu(par,ratio,sigma,vcoor,wN,Rot,reslist,kplot)
%--------------------------------------------------
%  df-sep-15: added output of corr coeff
%   df-oct-98  an economical version    df-may-98
%	visualize (plot) 1/(2R2R1-1) vs. sinTHETA^2
%--------------------------------------------------
if nargin<8, kplot=1; end		%default plot=on
Tx=par(1);
Dz2Dx=par(2);	
ee=Dz2Dx-1;
wT2=wN*wN*Tx*Tx;
nres=length(ratio);

costheta=(Rot*vcoor')';  		%direction cosines
sin2=1-costheta.^2;
theta=acos(costheta)*180/pi;
%------------calculate ratio ---------
calc=3/4/(1+wT2)*(1+ee*sin2.*wT2/(wT2+(1+ee/6)^2)./...
   (3+2*ee+(1+ee/3*(2-3*sin2)).^2).*...
   (4+3*ee+2/9*ee*ee-ee*sin2*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/...
   (wT2+(1+2*ee/3)^2))));
linear=3/4/(1+wT2)*(1+ee*sin2.*wT2/(wT2+1));
fitline=[[0;1],3/4/(1+wT2)*(1+ee*[0;1].*wT2/(wT2+1))];
thcrv=[0:2:180]';
sin2crv=sin(thcrv*pi/180).^2;
%------------calculate curve ---------
crv=3/4/(1+wT2)*(1+ee*sin2crv.*wT2/(wT2+(1+ee/6)^2)./...
   (3+2*ee+(1+ee/3*(2-3*sin2crv)).^2).*...
   (4+3*ee+2/9*ee*ee-ee*sin2crv*(1+(4+11/3*ee+19/18*ee*ee+5/54*ee*ee*ee)/...
   (wT2+(1+2*ee/3)^2))));
z=[sin2,ratio,sigma,calc,linear,theta,reslist];
thetacrv=[thcrv,crv];
mat2ascii('2D_data.txt',z)%record in file
mat2ascii('2D_data_calc.txt',thetacrv)
if kplot==0, return; end
%----------------plot-----------------




  figure(5)
  clf
  plot(sin2,ratio,'.b','MarkerSize',20);
  title('\fontsize{12}Axially Symmetric model');
  hold on
  for ii=1:nres,
    plot([sin2(ii) sin2(ii)],[ratio(ii)-sigma(ii) ratio(ii)+sigma(ii)]);
  end
  grid on 
  xlabel('\fontsize{12}sin^2(\theta)');
  ylabel('\fontsize{12}ratio');
  plot(sin2,calc,'b');
  plot(sin2,linear,'r');
  
  figure(6);
  clf
  plot(theta,ratio,'.b','MarkerSize',20);
  title('\fontsize{12}Axially Symmetric model');
  hold on
  for ii=1:nres,
    plot([theta(ii) theta(ii)],[ratio(ii)-sigma(ii) ratio(ii)+sigma(ii)]);
  end
  xlabel('\fontsize{12}\theta');
  ylabel('\fontsize{12}ratio');
  plot(thcrv,crv,'r');

  
  
  
  %------- plot difference between exp and calc ratio ----------------------
  
  p=[par(1) par(2)];
  theta_rad=theta*pi/180;
  ratio_calc=(calc_ratio_ax(p,theta_rad,wN))';   
  
diff=(ratio-ratio_calc)./sigma;  
co1=corrcoef(ratio,ratio_calc);  
corr=co1(1,2);
figure(7);
clf
subplot(2,1,1)
plot(ratio,ratio_calc,'.r','MarkerSize',20);
%plot(ratio,calc,'.r','MarkerSize',20);
title(['\fontsize{12}Axially Symmetric model (corr. coeff. = ',num2str(corr),')']);


sort_th=sort(ratio);
upper_bound=sort_th(nres);
lower_bound=sort_th(1);
hold on
grid on
xlabel('\fontsize{12}exp ratio');
ylabel('\fontsize{12}calc ratio');
axis([lower_bound*0.9 upper_bound*1.1 lower_bound*0.9 upper_bound*1.1])
plot(linspace(lower_bound*0.9,upper_bound*1.1,20),linspace(lower_bound*0.9,upper_bound*1.1,20));

sort_diff=sort(diff);
l_lim=sort_diff(1);
u_lim=sort_diff(nres);

figure(7)
subplot(2,1,2)
%JH change 8/06 for Matlab v7 %%%%%%%%%%%%%%%%
%bar(reslist,diff,'ro');
bar('v6',reslist,diff,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([min(reslist)-1 max(reslist)+1 l_lim u_lim])   %axis([0 reslist(nres) l_lim u_lim])
title('\fontsize{12}Axially Symmetric model');
grid on
xlabel('\fontsize{12}residue number');
ylabel('\fontsize{12}(exp ratio - calc ratio)/\sigma')

xxx=[reslist,diff];
save diff_ax.txt xxx -ascii
  
  
  %---Quality factor-----
num=(ratio-ratio_calc).^2;
num2=mean(num);
den=ratio-mean(ratio);
den1=den.^2;
den2=mean(den1);
den3=2*den2;

%den=ratio.^2;
%den1=mean(den);
%den2=2*den1;

Q=sqrt(num2/den3);





  
  
  
%----------------------------------------
return
%=====================================
