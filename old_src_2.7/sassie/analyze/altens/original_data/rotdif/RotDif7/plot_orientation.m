function [Q,ratio_th,corr]=plot_orientation(Dx,Dy,Dz,alpha,beta,gamma,freq,r2,r1,NH,reslist,r3,wN_full,l_lim,u_lim,flag)

%  df-sep-15: added output of corr coeff
% plot ratio versus alpha and beta for each NH vectors
% fully anisotropic case-
% alpha describes orientation of NH with respect to z axis
% beta describes orientation of NH in the xy plane
% ow-University of Maryland-2002


par0=[Dx Dy Dz];

%------ generate map of theta and phi with respect to rotation-diffusion axes -------------

thetamap=linspace(0,180,50);
thetamap=thetamap(:);
phimap=linspace(-180,180,50);
phimap=phimap(:);

l_theta=length(thetamap);
l_phi=length(phimap);

thetamap_rad=thetamap*pi/180.0;
phimap_rad=phimap*pi/180.0;





%---- map ratio, alpha and beta --------
rot_mat=rotation_matrix(alpha,beta,gamma);

for i=1:l_phi
    for j=1:l_theta
        coor_pdb(1,1)=sin(thetamap_rad(i))*cos(phimap_rad(j));
        coor_pdb(1,2)=sin(thetamap_rad(i))*sin(phimap_rad(j));
        coor_pdb(1,3)=cos(thetamap_rad(i));
        ratio_th1(i,j)=calc_ratio(rot_mat,coor_pdb,wN_full,par0);
        
    end
end



%---- determine alpha and beta for experimental data (vectors comming from pdb) ---------


[ratio,sigma,vcoor,rlist]=r2r1prep(r2,r1,NH,reslist,r3);
coord_r=(rot_mat*(vcoor'))'; %---- transform coordinate in the rot. diff frame -----

l_coor=length(coord_r(:,1));
nres=length(ratio);

for i=1:l_coor
    theta_exp(i)=acos(coord_r(i,3)/(sqrt(coord_r(i,1)^2+coord_r(i,2)^2+coord_r(i,3)^2)));
    phi_exp(i)=atan2(coord_r(i,2),coord_r(i,1));
end

theta_exp=theta_exp';
phi_exp=phi_exp';

x=theta_exp*180/pi;
y=phi_exp*180.0/pi;

comp_table=[x,y,ratio,sigma];






%----- plot difference between experimental and calculated data -------



ratio_th=calc_ratio(rot_mat,coord_r,wN_full,par0);
diff=(ratio-ratio_th)./sigma;
diff1=(ratio-ratio_th);
count=1:1:l_coor;




%---- Quality factor ------
num=(ratio-ratio_th).^2;
num2=mean(num);
den=ratio-mean(ratio);
den1=den.^2;
den2=mean(den1);
den3=2*den2;


%den=ratio.^2;
%den1=mean(den);
%den2=2*den1;


Q=sqrt(num2/den3);








%-----------------------------------------------
%----- plot the ratio versus exp ratio ----------
%-----------------------------------------------

figure(8)
clf
subplot(2,1,1)
hold on
grid on
co1=corrcoef(ratio,ratio_th);
corr=co1(1,2);
plot(ratio,ratio_th,'.r','MarkerSize',20);
title(['\fontsize{12}Fully Anisotropic model (corr. coeff. = ',num2str(corr),')']);


sort_th=sort(ratio);
upper_bound=sort_th(nres);
lower_bound=sort_th(1);

hold on
grid on
xlabel('\fontsize{12}exp ratio');
ylabel('\fontsize{12}calc ratio');
axis([lower_bound*0.9 upper_bound*1.1 lower_bound*0.9 upper_bound*1.1]);
plot(linspace(lower_bound*0.9,upper_bound*1.1,20),linspace(lower_bound*0.9,upper_bound*1.1,20));




%---- axis size ----

if(flag==2),
    figure(8)
    subplot(2,1,2)
    %JH change 8/06 for Matlab v7 %%%%%%%%%%%%%%%%
    %bar(rlist,diff,'r')
    bar('v6',rlist,diff,'r')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grid on;
    %axis([0 max(rlist) l_lim u_lim])
    axis([min(rlist)-1 max(rlist)+1 l_lim u_lim])
    title('\fontsize{12}Fully Anisotropic model');
    xlabel('\fontsize{12}residue number');
    ylabel('\fontsize{12}(exp ratio - calc ratio)/\sigma')
else
    l_lim=min(diff);
    u_lim=max(diff);
    figure(8)
    subplot(2,1,2)
    %bar(rlist,diff,'ro')
    bar('v6',rlist,diff,'r')
    grid on;
    %axis([0 max(rlist) l_lim u_lim])
    axis([min(rlist)-1 max(rlist)+1 l_lim u_lim])
    title('\fontsize{12}Fully Anisotropic model');
    xlabel('\fontsize{12}residue number');
    ylabel('\fontsize{12}(exp ratio - calc ratio)/\sigma')
end

xxx=[rlist,diff];
save diff_fan.txt xxx -ascii

%-----------------------------------------
%------ plot theoretical surface ---------
%-----------------------------------------

figure(9);
clf
p=surf(phimap,thetamap,ratio_th1,'FaceLighting','phong','EdgeColor','none');
set(p,'FaceAlpha',0.6);
grid on;
title('\fontsize{12}Fully Anisotropic model');
colormap(gray);
shading interp;
camlight;
lighting gouraud;


%hold on
%plot3(y,x,ratio,'.g','MarkerSize',20);
%hold on
  for ii=1:nres,
      
    if  (diff1(ii)>=0.0),
        hold on
        plot3([y(ii) y(ii)],[x(ii) x(ii)],[ratio(ii) ratio_th(ii)],'Color','r');
        hold on
        plot3(y(ii),x(ii),ratio(ii),'.r','MarkerSize',20);
    else
    hold on  
    plot3([y(ii) y(ii)],[x(ii) x(ii)],[ratio(ii) ratio_th(ii)],'Color','b');
    hold on
    plot3(y(ii),x(ii),ratio(ii),'.b','MarkerSize',20);
end
end

xlabel('\fontsize{12}\phi');
ylabel('\fontsize{12}\theta');
zlabel('\fontsize{12}ratio');


%-------- plot the projection of the 3D surface in 2D ----------
%-------- two extreme case : phi=0 and phi=90 ------------------


%-------- first case ------------------------

thetamap=linspace(0,180,50);

l_theta=length(thetamap);
thetamap_rad=thetamap*pi/180.0;

phi=0.0;
tab_theta=[];
for j=1:l_theta
        coor_pdb(1,1)=sin(thetamap_rad(j))*cos(phi);
        coor_pdb(1,2)=sin(thetamap_rad(j))*sin(phi);
        coor_pdb(1,3)=cos(thetamap_rad(j));
        ratio_th1=calc_ratio(rot_mat,coor_pdb,wN_full,par0);
        tab_theta=[tab_theta;ratio_th1];
end
figure(10)
clf
plot(thetamap(1,:),tab_theta(:,1),'--r','LineWidth',2);
title('\fontsize{12}Fully Anisotropic model');
xlabel('\fontsize{12}\theta');
ylabel('\fontsize{12}ratio');

%------ plot it on surface -------------
phi_tab=zeros(50,1).*0.0;
map_theta=thetamap';                %DF:stopped printout
figure(9)
hold on
plot3(phi_tab(:,1),map_theta(:,1),tab_theta(:,1),'--r','LineWidth',2);




%-------- second case ------------------------

phi=90.0*pi/180.0;
tab_theta2=[];
for j=1:l_theta
        coor_pdb(1,1)=sin(thetamap_rad(j))*cos(phi);
        coor_pdb(1,2)=sin(thetamap_rad(j))*sin(phi);
        coor_pdb(1,3)=cos(thetamap_rad(j));
        ratio_th2=calc_ratio(rot_mat,coor_pdb,wN_full,par0);
        tab_theta2=[tab_theta2;ratio_th2];
end


figure(10)
hold on
plot(thetamap(1,:),tab_theta2(:,1),'--g','LineWidth',2);

figure(10)
hold on
plot(comp_table(:,1),comp_table(:,3),'.b','MarkerSize',20);
for ii=1:nres,
    plot([comp_table(ii,1) comp_table(ii,1)],[comp_table(ii,3)-comp_table(ii,4) comp_table(ii,3)+comp_table(ii,4)]);
end

%--------- record -----------
tab3D_calc=[thetamap(1,:)',tab_theta(:,1),tab_theta2(:,1)];
tab3D=[comp_table(:,1),sigma,comp_table(:,3)];
mat2ascii('3D_data.txt',tab3D);%record in file
mat2ascii('3D_data_calc.txt',tab3D_calc);


%------ plot it on surface -------------
phi_tab=ones(50,1).*-90.0;
map_theta=thetamap';                 %DF:stopped printout  
figure(9)
hold on
plot3(phi_tab(:,1),map_theta(:,1),tab_theta2(:,1),'--g','LineWidth',2);

