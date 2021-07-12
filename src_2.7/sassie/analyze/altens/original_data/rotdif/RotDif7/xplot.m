function ratio=plot_orientation3(Dx,Dy,Dz,alpha,beta,gamma,freq,r2,r1,NH,reslist,r3)


wNfull=freq*0.1013756*2*pi*1.e06;
par0=[Dx Dy Dz];

%------generate map of theta and phi with respect to rot diff axes

thetamap=[-180:10:180];
%thetamap=thetamap(:);
phimap=[-180:10:180];
%phimap=phimap(:);

l_theta=length(thetamap);
l_phi=length(phimap);

thetamap_rad=thetamap*pi/180.0;
phimap_rad=phimap*pi/180.0;

%---- generate coordinates in the pdb frame --------
rot_mat=rotation_matrix(alpha*pi/180.0,beta*pi/180.0,gamma*pi/180.0);
k=1;
for i=1:l_theta
    for j=1:l_phi
        coor_pdb(1,1)=sin(thetamap_rad(i))*cos(phimap_rad(j));
        coor_pdb(1,2)=sin(thetamap_rad(i))*sin(phimap_rad(j));
        coor_pdb(1,3)=cos(thetamap_rad(i));
        ratio(i,j)=calc_ratio(rot_mat,coor_pdb,wNfull,par0);
    end
end



%--- plot th data in 3D ---------------

figure(1);
surf(thetamap,phimap,ratio,'EdgeColor','none');
camlight left; 
lighting gouraud
view(-15,65)

xlabel('phi');
ylabel('theta');
zlabel('ratio');

%view(3);axis tight; grid on;alpha(0.8);
%camlight headlight;
lighting gouraud;
hold off


%---- plot exp data in 3D ---------

[ratio,sigma,vcoor,rlist]=r2r1prep(r2,r1,NH,reslist,r3);
coord_r=(rot_mat*(vcoor'))';

l_coor=length(coord_r(:,1));

for i=1:l_coor
    theta_exp(i)=acos(coord_r(i,3)/(sqrt(coord_r(i,1)^2+coord_r(i,2)^2+coord_r(i,3)^2)));
    phi_exp(i)=atan(coord_r(i,2)/coord_r(i,1));
end

theta_exp=theta_exp';
phi_exp=phi_exp';
x=theta_exp*180/pi;
y=phi_exp*180.0/pi;

%plot3(x,y,ratio,'or','MarkerSize',8);
%plot3(y,x,ratio,'or','MarkerSize',8);
 
%stem3(y,x,ratio,'ok','fill');
%stem3(x,y,ratio,'ok','fill');


%----- plot diff between exp and theoretical data-------


ratio_th=calc_ratio(rot_mat,coord_r,wNfull,par0);

diff=abs(ratio-ratio_th);
diff2=(ratio-ratio_th);
count=1:1:l_coor;
length(rlist)
length(ratio)
nres=length(ratio);
figure(1)
hold on

plot3(y,x,ratio,'or','MarkerSize',4);
hold on
  for ii=1:nres,
    plot3([y(ii) y(ii)],[x(ii) x(ii)],[ratio_th(ii) ratio(ii)]);
end




figure(2);
plot(rlist,diff,'ro');

figure(4)
grid on
bar3(rlist,diff);




