function [thetamap,tab_theta,tab_theta2]=plot_orientation_demo(Dx,Dy,Dz,alpha,beta,gamma,freq,r2,r1,NH,reslist,r3,wN_full)

% plot ratio versus alpha and beta for each NH vectors
% fully anisotropic case-
% theta describes orientation of NH with respect to z axis
% phi describes orientation of NH in the xy plane
% ow-University of Maryland-2002


if nargin<13, wN_full=600.13*0.1013756*2*pi*1.e06;	end
if nargin<12, r3=[]; end
if nargin<11, reslist=[];	end		
if nargin<10, NH=[];	end		
if nargin<9, r1=[];	end
if nargin<8, r2=[];	end
if nargin<7, freq=600.13;	end
if nargin<6, gamma=42;	end
if nargin<5, beta=33;	end
if nargin<4, alpha=25;	end






par0=[Dx Dy Dz];

%------ generate map of theta and phi with respect to rotation-diffusion axes -------------



%---- map ratio, alpha and beta --------
rot_mat=rotation_matrix(alpha,beta,gamma);





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
thetamap=thetamap';
return


