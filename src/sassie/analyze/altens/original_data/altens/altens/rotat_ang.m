function rotat_ang(alpha1,beta1,gamma1,alpha2,beta2,gamma2,win_lin)


%---- if under windows, win_lin=0 ; if under linux win_lin=1

%-------- relaxation angles ------
alpha1_rad=alpha1*pi/180.0;
beta1_rad=beta1*pi/180.0;
gamma1_rad=gamma1*pi/180.0;


%-------- RDC angles -------------
alpha2_rad=[];beta2_rad=[];gamma2_rad=[]
alpha2_rad(1)=alpha2*pi/180.0;
beta2_rad(1)=beta2*pi/180.0;
gamma2_rad(1)=gamma2*pi/180.0;


%-------- explore all possibilities ---------
alpha2_rad(2)=alpha2_rad(1);
beta2_rad(2)=beta2_rad(1);
gamma2_rad(2)=gamma2_rad(1)+pi;

alpha2_rad(3)=pi+alpha2_rad(1);
beta2_rad(3)=pi-beta2_rad(1);
gamma2_rad(3)=pi-gamma2_rad(1);

alpha2_rad(4)=pi+alpha2_rad(1);
beta2_rad(4)=pi-beta2_rad(1);
gamma2_rad(4)=2*pi-gamma2_rad(1);








rot_mat1=rotation_matrix(alpha1_rad,beta1_rad,gamma1_rad);

tab_eul=[];
for ii=1:4
rot_mat2=rotation_matrix(alpha2_rad(ii),beta2_rad(ii),gamma2_rad(ii));


%res1=rot_mat1*[1;0;0]
%res2=rot_mat2*[1;0;0]

C=(rot_mat2*inv(rot_mat1));
%C=inv(rot_mat1)*rot_mat2


%res3=C*res1
e_ang=eulangMODIFIED(C);
tab_eul=[tab_eul;[e_ang(1),e_ang(2),e_ang(3)]];
th1=acos((dot(rot_mat1(:,1),rot_mat2(:,1))/(norm(rot_mat1(:,1))*norm(rot_mat2(:,1)))))*180/pi
th2=acos((dot(rot_mat1(:,2),rot_mat2(:,2))/(norm(rot_mat1(:,2))*norm(rot_mat2(:,2)))))*180/pi
th3=acos((dot(rot_mat1(:,3),rot_mat2(:,3))/(norm(rot_mat1(:,3))*norm(rot_mat2(:,3)))))*180/pi
end


fprintf('---------\n\n\n')
fprintf('alpha=%6.3f  beta=%6.3f  gamma=%6.3f\n\n\n',e_ang(1),e_ang(2),e_ang(3));




if (win_lin==1),
    %tab_solve=[];
    tab_lsq=[];
    for ii=1:4
     
   rot_mat2=rotation_matrix(alpha2_rad(ii),beta2_rad(ii),gamma2_rad(ii));
        C=rot_mat2*inv(rot_mat1);
    %C=inv(rot_mat1)*rot_mat2

%------- solving equations with least square procedure -------------

x0=[10*pi/180,2*pi/180.0,10*pi/180];
%options=optimset('Display','off','MaxFunEvals',50000,'MaxIter',10000,'TolFun',1e-15,'TolX',1e-15);

%x=fsolve(@myfun,x0,options,C);


%tab_solve=[tab_solve;[x(1,1)*180/pi,x(1,2)*180/pi,x(1,3)*180/pi]];

lb=[0 0 0];
ub=[2*pi 2*pi 2*pi];

options=optimset('Display','off','TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000,'LargeScale','on');

x_ax=lsqnonlin(@myfun2,x0,lb,ub,options,C);


rotation_matrix(x_ax(1,1),x_ax(1,2),x_ax(1,3));



tab_lsq=[tab_lsq;[x_ax(1,1)*180/pi,x_ax(1,2)*180/pi,x_ax(1,3)*180/pi]];
end
%tab_solve
tab_lsq
end
tab_eul
return