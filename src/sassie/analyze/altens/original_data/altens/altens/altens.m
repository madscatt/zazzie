function [table_all,Ds,R,corr,tens,S_tensor]=altens(D,NH,reslist,nsteps)
%----------------------------------------------------
%   Sept 2015: NOTE change in the order of input parameters!!!
%   modif: df-sep_15 df-apr_15  df-dec-11  df-jul-07
% df-sep-15: changed the order of input parameters (now can use all residues as default)
%            added to S_tensor to output
%            included euler2all selection of the angles from 0 to 180deg
%            renamed internal params to be consistent with their
%            meaning: tens --> eigvectors; S_non_diag --> S_tensor
%            renamed prim_angles -->prin_angles
% df-apr-15: added output tens and corr
%
% determine the alignment tensor extracted from 
% residual dipolar couplings
% derived from the previous version tens_align but now,
% we generate exactly the number of MC step and reject those which are not in the 
% boundary conditions
%   INPUT:
% D : n x 2 array; experimental residual dipolar couplings (first column=residue number;sec column:RDC)
% NH : n x 4 array of NH vector coordinates in the pdb frame [res#,x,y,z]
% reslist : list of residues to consider (default: all residues)
% nsteps : number of Monte Carlo steps for error report (default = 500)
%
% it can be usefull to use a combination of autopick_ellipsoid and res_dip to produce D
% currently the error in RDC is hard-wired (1 Hz?)
%
%DF: output: table_all = [res#, Dexp, vcoor, dir_cos]
%            Ds = [res#, Dexpt, Dcalc, Dexp-Dcalc]
%            R -- quality factor
%DF:         corr -- corr.coeff
%            tens=[Sx,Sy,Sz,alpha,beta,gamma];   
%            S_tensor -- the complete alignment tensor [3x3]
%---------------------------------------------------------


%=========================================================
%================= Initial calculation ===================
%=========================================================
format compact
if nargin < 4, nsteps = 500; end            %default: 500 MC steps
if nargin <3, reslist=[]; end               %default: take all residues
if isempty(reslist), reslist=D(:,1); end    %default: take all residues

%--------- prepare input --------------


[Di,vcoor,rlist]=prep2(D,NH,reslist);

max_di=max(Di);
min_di=min(Di);

fprintf('minimum dipolar coupling : %7.3f  ;  maximum dipolar coupling : %7.3f \n\n',min_di,max_di);
 
table_all=[rlist,Di,vcoor];

%-------- plot residual dipolar coupling for each residues ---------------
figure(1)
clf
bar(table_all(:,1),table_all(:,2))
axis([0 max(table_all(:,1)) min(Di)-1 max(Di)+1])
title('\fontsize{14}RDCs vs residue number')
xlabel('\fontsize{14}Residue number')
ylabel('\fontsize{14}RDC')

nres=size(table_all,1);
disp([num2str(nres),' residues selected']);
fprintf('\n');

sortrows(table_all,2);
fprintf('==============================================================================\n');
fprintf('===================== Alignment Tensor characteritics =========================\n');
fprintf('==============================================================================\n\n');

%--------- calculate direction cosine --------------

z=dir_cos(vcoor);

%--------- add it to the table -----------


table_all=[rlist,Di,vcoor,z];

%-------- generate matrix A --------------

[mat_A,mat_A_inv,t,v,con_num]=gen_A(table_all(:,6:8));


%-------- calculate matrix x -------------

[eigvectors,S,S_tensor]=x_tensor3(mat_A_inv,table_all(:,2))
%[tens,S,S_non_diag]=x_tensor3(mat_A_inv,table_all(:,2))
%fprintf('Syy=%6.3f Szz=%6.3f  Sxy=%6.3f  Sxz=%6.3f  Syz=%6.3f\n\n\n',S_non_diag(2,2),S_non_diag(3,3),S_non_diag(1,2),S_non_diag(1,3),S_non_diag(2,3));

fprintf('\n');
fprintf('Magnitude\n');
fprintf('---------\n\n');
fprintf('Sxx''=%6.3f  Syy''=%6.3f  Szz''=%6.3f\n\n',S(1),S(2),S(3));



%-------- extract euler angles ------------

e_ang=euler_angle(eigvectors');
[all_eulers, select_eulers] = euler2all(e_ang(1),e_ang(2),e_ang(3));   %df-sep-15
e_ang=select_eulers;                                                   %df-sep-15     

fprintf('Orientation\n');
fprintf('---------\n\n')
fprintf('alpha=%6.3f  beta=%6.3f  gamma=%6.3f\n\n',e_ang(1),e_ang(2),e_ang(3));


prin_angle=e_ang;
prin_param=S;

tens=[S,e_ang];     %DF: apr 2015 

%============ test the quality of the fit ==========
%------------ extract Sij and recalculate Dij -----------

[Syy,Szz,Sxy,Sxz,Syz]=extract_par(S_tensor);
new_D=recalc_D(Syy,Szz,Sxy,Sxz,Syz,mat_A);

N=[new_D,Di];
Ds=[table_all(:,1),Di,new_D,Di-new_D];      %DF:report the exp and calc values and their difference
figure(2)
clf
stem(Ds(:,1),Ds(:,2)-Ds(:,3),'fill','r');   %DF: corrected the sign of the difference
%stem(table_all(:,1),N(:,1)-N(:,2),'fill','r');
grid on
title('\fontsize{14}Difference between experimental & recalculated RDC (Hz)');
xlabel('\fontsize{14}Residue number');
ylabel('\fontsize{14}RDCexp - RDCcalc (Hz)');


%---------- calculate R factor -------------

corr = corrcoef(Di,new_D);
corr = corr(1,2);                   %DF: apr 2015
R=calc_R(Di,new_D);
An_param=calc_ani(prin_param);
fprintf('corr.coef = %5.3f ; R-factor = %6.3f ;  asymmetry parameter = %6.3f \n\n',corr,R,An_param);

%---------- Plot correlation -------------

figure(3)
clf
plot_diff(new_D,Di,R);        %DF: the program plot_diff inverts the order, such that Di is abscissa 

%return

%===============================================================
%======================= Monte Carlo ===========================
%===============================================================

res=input('Monte Carlo error? ( [1]---> yes; 0 --->no) ==> ');
if isempty(res), res=1; end
if (res==1)
    
    fprintf('\n');
    fprintf('calculation in progress........\n\n');
    
    
    %=================== Monte Carlo =================
    e_set=NaN*ones(nsteps,3);
    S_set=NaN*ones(nsteps,3);
    e_ang_full=[];
    tab_param=[];
    tab_ang=[];
    n_gen=0;
    k=0;
    while (k~=nsteps),
        n_gen=n_gen+1;
        table_new=mc_table(Di);



    %-------- calculate matrix x -------------

    [eigvalues_mc,S_mc,S_tensor_mc]=x_tensor3(mat_A_inv,table_new(:,1));
    
  
    %-------- extract euler angles ------------

    e_ang=euler_angle(eigvalues_mc');
    [all_eulers, select_eulers] = euler2all(e_ang(1),e_ang(2),e_ang(3));   %df-sep-15
    e_ang=select_eulers;                                                   %df-sep-15 
    
    e_ang_full=[e_ang_full;e_ang];
    
    %------- check if S principal parameters are meaningful------------
    
        if ((S_mc(1)<=(prin_param(1)+5.0))&(S_mc(1)>=(prin_param(1)-5.0))&(S_mc(2)<=(prin_param(2)+5.0))&(S_mc(2)>=(prin_param(2)-5.0))&...
                (S_mc(3)<=(prin_param(3)+5.0))&(S_mc(3)>=(prin_param(3)-5.0))&(e_ang(1)<=(prin_angle(1)+60.0))&(e_ang(1)>=(prin_angle(1)-60.0))...
                &(e_ang(2)<=(prin_angle(2)+60.0))&(e_ang(2)>=(prin_angle(2)-60.0))&(e_ang(3)<=(prin_angle(3)+60.0))&(e_ang(3)>=(prin_angle(3)-60.0))),
            k=k+1;
         
            tab_param=[tab_param;[S_mc(1),S_mc(2),S_mc(3)]];
            tab_ang=[tab_ang;[e_ang(1),e_ang(2),e_ang(3)]];
            test_progress(k,nsteps);
            
        end
end
    
    
 
%------- useful prameters -------------------

use_param=size(tab_param,1);

%------- useful angles ----------------------

use_angles=size(tab_ang,1);
   
   
%------- calculate mean and stdv for principal parameters -------------------- 
Sxx=mean(tab_param(:,1));
std_Sxx=std(tab_param(:,1));

Syy=mean(tab_param(:,2));
std_Syy=std(tab_param(:,2));

Szz=mean(tab_param(:,3));
std_Szz=std(tab_param(:,3));
   
   
%------- calculate mean and stdv for angles -----------------------------------    
res_alpha=mean(tab_ang(:,1));
sd_res_alpha=std(tab_ang(:,1));

res_beta=mean(tab_ang(:,2));
sd_res_beta=std(tab_ang(:,2));

res_gamma=mean(tab_ang(:,3));
sd_res_gamma=std(tab_ang(:,3));

%-------- plot distribution of parameters (comming from MC) ------------------

plot_mc(tab_ang,tab_param)

    
param=[Sxx Syy Szz];    
An_param=calc_ani(param);
error_ani=An_param*sqrt(((std_Sxx^2+std_Syy^2)/(Sxx-Syy)^2)+(std_Szz^2/Szz^2));

fprintf('\n');
fprintf('======================================================\n');
fprintf('========== Monte Carlo error Report ==================\n');
fprintf('======================================================\n\n');




fprintf('number of generated parameters=%6.0d \n\n',n_gen);

fprintf('Sxx=%6.3f+/-%6.3f  Syy=%6.3f+/-%6.3f  Szz=%6.3f+/-%6.3f\n\n',Sxx,std_Sxx,Syy,std_Syy,Szz,std_Szz);

fprintf('alpha=%6.3f+/-%6.3f  beta=%6.3f+/-%6.3f  gamma=%6.3f+/-%6.3f\n\n',res_alpha,sd_res_alpha,res_beta,sd_res_beta,res_gamma,sd_res_gamma);

fprintf('asymmetry parameter = %6.3f +/- %6.3f \n\n',An_param,error_ani);

end
