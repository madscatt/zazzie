function [tab_gen]=RotDif(freq,r1,r2,r3,NH,reslist,output,pdb_file,alg,mc_max,niter,guessMC_ax,guessMC_full)
%------------------------------------------------------------------------------------------------------------------------
%   df-sep-15   redesigned output and order of computations: isotropic model comes first and without asking)
%               FINAL reports contain real minima and MC-derived standard deviations
%               corrected computation of errors in rhombicity, now compute 
%               tauC, ani, and rhomb errors in one function
%               the program now reports corr coeff for axi and full models (modif r2r1visu and plot_orient..)
%               next thing to do: get both Q and corr coeffs included in the final report 
%               
%   ow-2002-University of Maryland
%   based on the program r2r1fit from david fushman
%
%
%   Optimization of isotropic, axially symmetric and fully anisotropic model
%   
%   CALL: RotDif(freq,r1,r2,r3,NH,reslist,output,pdb_file,alg,mc_max,niter,guessMC_ax,guessMC_full)
%
%   INPUT
%   ------
%   freq : 1H spectrometer frequency in MHz
%   r1,r2,r3 : R1, R2, NOE in 1/ms
%   reslist : list of residues to be taken for the analysis ([]=take all)
%   output : name of output file (input a string variable, e.g. 'tmp') --> the data will 
%            be written to file tmp_date_time.txt
%   alg : algorithm to use for least square search; 1=Levenberg-Marquardt (default); 0=newton
%   pdb_file: Name of pdb file in order to extract NH vectors. If not present, the user 
%               must provide NH vectors (vNH). If a pdb file is available and NH vector 
%               coordinates are to be taken from this file, put vNH=[];
%   guessMC_ax : 1*8 array containing boundary conditions to simulate starting 
%                 guess prior to data fitting (axially symmetric case) and boundary conditions
%                  (used in the fitting procedure for axially symmetric case) 
%   [lower bound for TAUx,upper bound for TAUX,lower bound for D/D,.....alpha,....beta]
%   
%   guessMC_full: same as guessMC_ax for Dxx, Dyy, Dzz, alpha, beta, gamma ----> 1*12 array
%	      
%   mc_max : number of MC step for error calculation (default=1000)
%   niter : number of simulated starting guess prior to optimization (default=50)
%
%
%   OUTPUT
%   ------
%
%   tab_gen: array containing [residue1, sigma1, exp ratio1, calc ratio1 for F.An. model, calc ratio1 for Ax.Sym model]
%                             [residue2, sigma2, exp ratio2, calc ratio2 for F.An. model, calc ratio2 for Ax.Sym model]  
%                                   |       |                   |                           |                   |
%                                   |       |                   |                           |                   |
%                                   |       |                   |                           |                   |
%                                   |       |                   |                           |                   |
%
%-------------------------------------------------------------------------------------------------------------------------
% If the three dynamic model (isotropic, axially symmetric, fully anisotropic) are considered
% 5 text files are saved in the current directory (to transfer data in origin)
% Each text file corresponds to a given figure
%
% 3D_data.txt : to plot figure 10 (exp points)
%	[theta; sigma; ratio] 
%
% 3D_data_calc.txt : to plot figure 10 (theoretical curves)
%	[theta calculated; red curve(upper values); green curve(lower values)]
%
% 2D_data.txt : to plot figures 5 and 6
%	[sin2,ratio; sigma; calculated ratio; ratio if linear; theta; reslist]
%
% 2D_data_calc.txt : to plottheoretical curve of figure 6
%	[theta; ratio calculated]
%
% tab_gen.txt : to plot figures 7 and 8
%	[reslist; sigma; ratio; ratio calculated with fully anisotropic model; ratio calculated with ax. sym. model]
%
%
%
%
%
%--------------------------------------------------------------------------------------------------------------------------


format compact
%--------check input, set some params-------


if nargin<13, guessMC_full=[0.5 12.0 0.5 12.0 0.5 12.0 0 pi 0.0 pi 0.0 pi];	end
if nargin<12,  guessMC_ax=[1.0 15.0 1.0 2.5 0 pi 0 pi]; end    %for oblate: [2.0 12.0 0.05 1.0 0 pi 0 pi];
if nargin<11, niter=50;	end		%default: number of initial values for MC estimation
if nargin<10, mc_max=1000;	end		%default: number of MC steps for statistical errors 
if nargin<9, alg=1;	end
if ((nargin<8)&(isempty(NH)))
    error('ERROR, unable to find at least one NH coordinate file')
end
from_pdb=0;
if nargin < 7, output = 'junk';  end        %default    %DF
if nargin < 6, reslist = [];  end        %default       %DF
if isempty(reslist), reslist=r1(:,1)'; end  %take all residues  %DF

if (isempty(NH))
    from_pdb=1;
    NH=pdb2nh(pdb_file,[],1);
end

wN=freq*(2710.5/26750)*2*pi/1000;
wNfull=freq*0.1013756*2*pi*1.e06;
flag_case=0;     %check for axially symmetric or fully anisotropic case
flag_visu=0;
fan_calc=[];ax_calc=[];

bound_full=guessMC_full;
bound_ax=guessMC_ax;

%----- save data -------
cl=clock;
fileout=strcat(output,'_',num2str(cl(2)),'_',num2str(cl(3)),'_',num2str(cl(4)),'_',num2str(cl(5)),'.txt');
fid_out=fopen(fileout,'w');


fprintf('\n\n Your Results will be saved in %s \n\n',fileout);


%----- record pdb filename and date --------
fprintf(fid_out,'date : %2.0d/%2.0d/%4.0d at %2.0d:%2.0d \n',cl(2),cl(3),cl(1),cl(4),cl(5));
if (from_pdb==1)
    fprintf(fid_out,'NH vectors extracted from : %s \n\n\n',pdb_file);
else
    fprintf(fid_out,'Read NH vectors provided by user \n\n\n');
end



%------------check output-------------------
if nargout>=4, kmap=1; else kmap=0; end
sdzres=NaN*ones(1,4);
%----------prepare working set of data-----------
[ratio,sigma,vcoor,rlist]=r2r1prep(r2,r1,NH,reslist,r3);
nres=length(ratio);
dfiso=nres-1;     %degree of freedom isotropic case
df=nres-4;        %degree of freedom axially symmetric
df_full=nres-6;   %degree of freedom fully anisotropic
fprintf('\n');
disp([num2str(nres),' residues selected']);
fprintf('\n')
init_tab=[rlist,sigma,ratio];




%------ extract and record unused residues -----
uvec=sort_unused(reslist,r1(:,1));
fmt=' %3.0d';
[nrows, ncols]  = size(uvec);
fprintf('Excluded residues : \n');
fprintf([repmat(fmt,1,ncols)],uvec');
fprintf(' \n');
fprintf(fid_out,'Excluded residues : \n\n');
%fprintf(fid_out,[repmat(fmt,1,ncols) '\n'], uvec');
fprintf(fid_out,[repmat(fmt,1,ncols)], uvec');      %df modif
fprintf(fid_out,'\n\n');




%------ save in output file --------
fprintf(fid_out,'%3.0d residues selected\n',nres);




%---- plot ratio vs residue number ------
figure(1);
clf
plot(rlist,ratio,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                                'MarkerSize',5)
grid on
title('\fontsize{12}ratio vs residue number');
xlabel('\fontsize{12}residue number');
ylabel('\fontsize{12}ratio');



%--------- show space sampling by NH vectors----------
figure(2);
clf
arrow3(zeros(size(vcoor,1),3),vcoor,'b',0.05,0.1);
grid on
light('Position',[0,-1,1]);
lighting gouraud;
axis square
title('\fontsize{12}Distribution of NH vectors');
%figure(2);
%clf
%draw_vect(vcoor);

%=======================================================================================
%                                                                                      =
% ===============================  ISOTROPIC MODEL =====================================
%                                                                                      =
%=======================================================================================

disp(' ');

tauc_app=calc_each(wN,ratio);
me_tauc_eff=mean(tauc_app);
sd_tauc_eff=std(tauc_app);
fprintf('mean of all apparent TAUc values: %6.4f \n',me_tauc_eff);
fprintf('std dev of all apparent TAUc values: %6.4f \n',sd_tauc_eff);

%-------------------find min for isotropic model -----------------
optiso=[0 1.e-8 1.e-8 0 0 0 0 0 0 0 0 0 0 10000];
[pariso,chi2iso]=r2r1fminiso(me_tauc_eff/2,me_tauc_eff*2,optiso,ratio,sigma,wN);

%------------ quality factor ----------------
Q=quality_iso(pariso,ratio,sigma,wN);

%-------------prepare iso-report-------------

d_iso=tauc2d_iso(pariso);

fprintf('\n');
fprintf('================= ISOTROPIC MODEL =====================\n');
fprintf('TAUc = %5.2f ns  Diso = %5.2f  ; Chi2 = %9.4f  ; Chi2/df = %8.4f \n',pariso,d_iso,chi2iso,chi2iso/dfiso);
fprintf('Quality Factor = %5.3f \n',Q);
fprintf('=======================================================\n');


%----- write to output file -------
fprintf(fid_out,'\n');
fprintf(fid_out,'====================== ISOTROPIC MODEL ================\n');
fprintf(fid_out,'TAUc = %5.2f ns Diso = %5.2f ; Chi2 = %9.4f  ; Chi2/df = %8.4f \n',pariso,d_iso,chi2iso,chi2iso/dfiso);
fprintf(fid_out,'Quality Factor = %5.3f \n',Q);
fprintf(fid_out,'=======================================================\n');


invcov_iso=invcovar_iso(d_iso,sigma,freq);
[me_mc,sd_mc]=r2r1conf_iso(d_iso,chi2iso/dfiso,invcov_iso,mc_max);
    
tau_err=(pariso/(me_mc(1)))*sd_mc(1);
d_err=(d_iso/(me_mc(1)*1e-7))*sd_mc(1)*1e-7;
   
   
fprintf('\n');
fprintf('*********** MONTE CARLO error report for Isotropic model ************* \n');
fprintf('TAUc =%5.2f +/- %5.2f ns ; Diso =%5.2f +/- %5.2f \n',pariso,tau_err,me_mc(1)*1e-7,sd_mc(1)*1e-7);
fprintf('********************************************************************** \n\n');
   
   
fprintf(fid_out,'\n');
fprintf(fid_out,'********* MONTE CARLO error report for Isotropic model ********** \n');
fprintf(fid_out,'TAUc =%5.2f +/- %5.2f ns ; Diso =%5.2f +/- %5.2f \n\n',pariso,tau_err,me_mc(1)*1e-7,sd_mc(1)*1e-7);
fprintf(fid_out,'***************************************************************** \n\n');
      

%--------- sampling tensor statistics-----------------
SS=s_tens(vcoor,1);
SSeig=eig(SS)';
ff=ss2ff(SS)';

disp('Distribution of NH vectors :');
disp(' ');
disp(['Sampling param KSI: ',num2str(ff2ksi(ff)),...
      '; tensor components: ',num2str(SSeig)]);
fprintf('\n');

%-------- sampling statistics on output --------------
fprintf(fid_out,'Distribution of vectors :');
fprintf(fid_out,'\n\n');
fprintf(fid_out,'Sampling param KSI:%8.5f ; Tensor components : %8.5f %8.5f %8.5f \n\n\n',ff2ksi(ff),SSeig(1),SSeig(2),SSeig(3));



%=======================================================================================
%                                                                                      =
% =======================  AXIALLY SYMMETRIC MODEL =====================================
%                                                                                      =
%=======================================================================================

res=input('axially symmetric (+isotropic) model? [1]-->yes; [0]-->no : ');
if isempty(res), res=1; end
if (res==1)
    flag_case=flag_case+1;

%---------------- Monte carlo distribution to create a set of initial values --------------

 
count=[1:1:niter];  %-- count number of synthetic data, used for drawing
trialrad=ones(niter,4);

%------ Array of random initial guess to search the best solution-----
for i=1:niter
    trialrad(:,1)=guessMC_ax(1)+(guessMC_ax(2)-guessMC_ax(1))*rand(niter,1);
    trialrad(:,2)=guessMC_ax(3)+(guessMC_ax(4)-guessMC_ax(3))*rand(niter,1);
    trialrad(:,3)=guessMC_ax(5)+(guessMC_ax(6)-guessMC_ax(5))*rand(niter,1);
    trialrad(:,4)=guessMC_ax(7)+(guessMC_ax(8)-guessMC_ax(7))*rand(niter,1);
end

fprintf('Search Range\n');
fprintf('%5.1f < TAUx  < %5.1f\n',guessMC_ax(1),guessMC_ax(2));
fprintf('%5.1f < D/D   < %5.1f\n',guessMC_ax(3),guessMC_ax(4));
fprintf('%5.1f < alpha < %5.1f\n',guessMC_ax(5)*180.0/pi,guessMC_ax(6)*180.0/pi);
fprintf('%5.1f < beta  < %5.1f\n',guessMC_ax(7)*180.0/pi,guessMC_ax(8)*180.0/pi);
fprintf('\n');


%------------ angles in degree-----------------

for i=1:niter
    trialdeg(:,1)=trialrad(:,1);
    trialdeg(:,2)=trialrad(:,2);
    trialdeg(:,3)=trialrad(:,3)*180.0/pi;
    trialdeg(:,4)=trialrad(:,4)*180.0/pi;
end





%----------- calculate difference between calculated and exp data -----

reslsq=ones(niter,5);  %for LVM algorithm


%----- bounding conditions in the fitting procedure-------

lb=[bound_ax(1) bound_ax(3) bound_ax(5) bound_ax(7)]; %lower bound for TAUc D/D alpha beta
ub=[bound_ax(2) bound_ax(4) bound_ax(6) bound_ax(8)];   %upper bound for TAUc D/D alpha beta


%----- Fitting procedure ---------------------------------
if (alg==1),    
    %options=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000,'LargeScale','on','LevenbergMarquardt','on');
    options=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',100000);
else
    options=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',100000);
end

fprintf('\n');  
fprintf('********************************************************************************************\n');
fprintf('****************  Starting Monte Carlo simulation with %.0f initial values  ****************\n',niter);
fprintf('****************               for axially symmetric model                  ****************\n');
fprintf('********************************************************************************************\n\n');
tic

fprintf('Calculation in progress\n');
fprintf('0%%');
for i=1:niter-6
    fprintf(' ');
end
fprintf('100%% \n');
    

%----begin the MC search-------
for i=1:niter
    x_ax=lsqnonlin(@rho2,trialrad(i,:),lb,ub,options,ratio,sigma,vcoor,wN); %---- LVM algorithm----- 
    errlsq_ax=rho(ratio,sigma,vcoor,wN,x_ax);          % Chi^2 for LVM algorithm
    
    
    %------ calc progress ---------
   
    fprintf('>');
         
    %---- array of results with Levenberg Marquardt-----
    reslsq(i,1)=x_ax(1);
    reslsq(i,2)=x_ax(2);
    reslsq(i,3)=x_ax(3);
    reslsq(i,4)=x_ax(4);
    reslsq(i,5)=errlsq_ax;
end
%reslsq
%pause
fprintf('\n');
%---- plot the target function for the different initial values------
figure(3);
clf
subplot(2,1,1);
hist(reslsq(:,5),20)
title('\fontsize{12}Axially Symmetric model');
xlabel('\fontsize{12}\chi^{2}');
ylabel('\fontsize{12}N','Rotation',0);

toc


%---reorder matrix of results to pick the best target function ---------
%---- in radians ----
sortlsqrad=sortrows(reslsq,5);


%---- calculate TAUc for the best fit --------
Tauclsq=tx2tc(sortlsqrad(1,1),sortlsqrad(1,2));

%------ prepare report --------
zres_ax=[sortlsqrad(1,:)];

%----- calculate Dxx, Dyy, Dzz -------
endDx=calcDx(zres_ax(1));
endDz=calc_Dz(zres_ax(2),zres_ax(1));

fprintf('\n');
fprintf('================ AXIALLY SYMMETRIC MODEL ======================== \n');
fprintf('Axially symmetric model, Final results : \n');
fprintf('TAUx = %5.2f ns ;  D/D = %5.2f  ; alpha = %4.0f deg  ;  beta = %4.0f deg\n',zres_ax(1),zres_ax(2),zres_ax(3)*180.0/pi,...
    zres_ax(4)*180./pi);
fprintf('TAUC = %5.2f ns ; Chi2 = %9.4f ; Chi2/df = %8.4f\n',Tauclsq,zres_ax(5),zres_ax(5)/df);
fprintf('Dxx=Dyy = %5.2f*10^7  Dzz = %5.2f*10^7 \n',endDx,endDz);
fprintf('================================================================== \n');



%------- save in output file -------------------
fprintf(fid_out,'\n\n');
fprintf(fid_out,'================== AXIALLY SYMMETRIC MODEL ====================== \n\n');
fprintf(fid_out,'Axially symmetric model, Final results : \n\n');
fprintf(fid_out,'TAUx = %5.2f ns ;  D/D = %6.2f  ; alpha = %4.0f deg  ;  beta = %4.0f deg\n \n',zres_ax(1),zres_ax(2),zres_ax(3)*180.0/pi,...
    zres_ax(4)*180./pi);
fprintf(fid_out,'TAUC = %5.2f ns ; Chi2 = %9.4f ; Chi2/df = %8.4f\n \n',Tauclsq,zres_ax(5),zres_ax(5)/df);
fprintf(fid_out,'Dxx=Dyy = %5.2f*10^7  Dzz = %5.2f*10^7 \n\n',endDx,endDz);
fprintf(fid_out,'================================================================== \n\n');


%------compare with the isotropic model------
%--------------F-test aniso vs. iso----------
disp(' ');
[P,F]=ftestcmp(chi2iso/dfiso,dfiso,zres_ax(5)/df,df);
fprintf('F-Statistics ANISO vs. ISO: F= %8.5f, P= %8.5e \n\n',F,P);


%----- save to output file -------
fprintf(fid_out,'F-Statistics ANISO vs. ISO: F= %8.5f, P= %8.5e \n\n',F,P);


%--------------draw diff. axis---------------
figure(2);
alr=zres_ax(3);ber=zres_ax(4);
Daxis=2*[sin(ber)*cos(alr),sin(ber)*sin(alr),cos(ber)];
hold on
draw_vect([Daxis;-Daxis],'r',2.5);
axis square



%--------------standard errors in TAU & D/D using inverse covariance matrix procedure -------------
kerr=input('calculate errors (axially symmetric model)? ([1]--->yes,0--->no) ==> ');
if isempty(kerr), kerr=1; end

if kerr==1,
   invCov=r2r1covMC2(zres_ax(1:4),ratio,sigma,vcoor,wN);
   [me_mc,sd_mc]=r2r1mc(zres_ax(1:4),zres_ax(5)/df,invCov,mc_max,1);
   tc_mc=tx2tc(me_mc(1),me_mc(2));      
   sd_tc_mc=tc_mc*sqrt((sd_mc(1)/me_mc(1))^2+(sd_mc(2)/(me_mc(2)+2))^2);    % error in Tauc
   sdzres(1:4)=sd_mc;
   
   endDx2=calcDx(me_mc(1));
   endDz2=(calc_Dz(me_mc(2),me_mc(1)));
   errDx=(endDx2/me_mc(1))*sd_mc(1);                                %error in Dxx, Dyy
   errDz=endDz2*sqrt((sd_mc(2)/me_mc(2))^2+(sd_mc(1)/me_mc(1))^2);         %error in Dz
   
   fprintf('\n');
   fprintf('********* MONTE CARLO error report for Axially Symmetric Model ************* \n');
   fprintf('TAUx =%5.2f +/-%5.2f ns ; D/D =%5.2f +/-%5.2f ;\n',me_mc(1),sd_mc(1),me_mc(2),sd_mc(2));   
   fprintf('alpha =%4.0f +/-%4.0f deg ;  beta =%4.0f +/-%4.0f deg \n',me_mc(3)*180.0/pi,sd_mc(3)*180.0/pi,me_mc(4)*180.0/pi,sd_mc(4)*180.0/pi);
   fprintf('**************************************************************************** \n');
   
   fprintf('\n');
   fprintf('********* FINAL report for Axially Symmetric Model ************* \n');
   fprintf('Dxx=Dyy = %5.2f +/- %5.2f*10^7   Dz = %5.2f +/- %5.2f*10^7 \n',endDx,errDx,endDz,errDz);
   fprintf('alpha =%4.0f +/-%4.0f deg ;  beta =%4.0f +/-%4.0f deg \n',zres_ax(3)*180.0/pi,sd_mc(3)*180.0/pi,zres_ax(4)*180.0/pi,sd_mc(4)*180.0/pi);
   fprintf('TAUc =%5.2f +/-%5.2f ns ; D/D =%5.2f +/-%5.2f ; TAUx =%5.2f +/-%5.2f ns \n',Tauclsq,sd_tc_mc,zres_ax(2),sd_mc(2),zres_ax(1),sd_mc(1));   
   fprintf('**************************************************************************** \n\n');   
   
   
   fprintf(fid_out,'\n');
   fprintf(fid_out,'********* MONTE CARLO error report for Axially Symmetric Model ************* \n');
   fprintf(fid_out,'TAUx =%5.2f +/-%5.2f ns ; D/D =%5.2f +/-%5.2f ;\n',me_mc(1),sd_mc(1),me_mc(2),sd_mc(2));   
   fprintf(fid_out,'alpha =%4.0f +/-%4.0f deg ; beta =%4.0f +/-%4.0f deg \n',me_mc(3)*180.0/pi,sd_mc(3)*180.0/pi,me_mc(4)*180.0/pi,sd_mc(4)*180.0/pi);
   fprintf(fid_out,'**************************************************************************** \n\n');

   fprintf(fid_out,'\n');
   fprintf(fid_out,'********* FINAL report for Axially Symmetric Model ************* \n');
   fprintf(fid_out,'Dxx=Dyy = %5.2f +/- %5.2f*10^7   Dz = %5.2f +/- %5.2f*10^7 \n',endDx,errDx,endDz,errDz);
   fprintf(fid_out,'alpha =%4.0f +/-%4.0f deg ; beta =%4.0f +/-%4.0f deg \n',zres_ax(3)*180.0/pi,sd_mc(3)*180.0/pi,zres_ax(4)*180.0/pi,sd_mc(4)*180.0/pi);
   fprintf(fid_out,'TAUc =%5.2f +/-%5.2f ns ; D/D =%5.2f +/-%5.2f ; TAUx =%5.2f +/-%5.2f ns \n',Tauclsq,sd_tc_mc,zres_ax(2),sd_mc(2),zres_ax(1),sd_mc(1));   
   fprintf(fid_out,'**************************************************************************** \n');   

   
end



%----- plot target function for best results------

res3d=input('plot Chi^2 versus alpha and beta? ( [1]---> yes; [0]--->no) ==> ');
if isempty(res3d), res3d=1; end
if (res3d==1)
param=ones(1,4);
alphamap=linspace(0,guessMC_ax(6)*180.0/pi,50);
betamap=linspace(0,guessMC_ax(8)*180.0/pi,50);

la=length(alphamap);
lb=length(betamap);

for i=1:la
   for j=1:lb
   param=[zres_ax(1),zres_ax(2),alphamap(i)*pi/180.0,betamap(j)*pi/180.0];
   tf(i,j)=rho(ratio,sigma,vcoor,wN,param);
end
end


figure(4);
clf
surf(betamap,alphamap,tf,'FaceLighting','phong','EdgeColor','none');
colorbar;
grid on;
alpha(0.8);
shading interp
camlight;
lighting gouraud;
hold on;
contour(betamap,alphamap,tf,60)
%pcolor(betamap,alphamap,tf);
%shading interp;
hold on;
stem3(zres_ax(4)*180.0/pi,zres_ax(3)*180.0/pi,zres_ax(5),'ro','fill');

if (guessMC_ax(6)~=guessMC_ax(8))
    if ((zres_ax(4)*180.0/pi)>180),
        hold on
        stem3((zres_ax(4)*180.0/pi)-180,(zres_ax(3)*180.0/pi)+180,zres_ax(5),'ro','fill');
    end
    if ((zres_ax(3)*180.0/pi)>180),
        hold on
        stem3(180-(zres_ax(4)*180.0/pi),180-(zres_ax(3)*180.0/pi),zres_ax(5),'ro','fill');
    end
    hold on
    stem3(180-(zres_ax(3)*180.0/pi),(zres_ax(4)*180.0/pi)+180,zres_ax(5),'ro','fill');
end
%legend('\fontsize{12}\chi^{2} minimum',1);
xlabel('\fontsize{12}\beta');
ylabel('\fontsize{12}\alpha');
zlabel('\fontsize{12}\chi^{2}','Rotation',0);
title('\fontsize{12}\chi^{2} versus \alpha and \beta (ax.sym.model)');


end


%--------------------visualize points?-----------------
kvis=input('visualize? ([1]-yes,0-no) ==> ');
if isempty(kvis), kvis=1; end
if kvis==1, 
flag_visu=flag_visu+1;
zvis=[];
fitline=[];
thetacrv=[];
par=[zres_ax(1) zres_ax(2) zres_ax(3) zres_ax(4)];
rotagrid=[sin(zres_ax(4))*cos(zres_ax(3)) sin(zres_ax(4))*sin(zres_ax(3)) cos(zres_ax(4))];

  
  [zvis,fitline,thetacrv,l_lim,u_lim,Q_ax,corr_ax,ax_calc]=r2r1visu(par,ratio,sigma,vcoor,wN,...
      rotagrid,rlist);
  
fprintf('\n');  
fprintf('Quality Factor for Axially Symmetric Model = %5.3f ; Corr.coeff = %7.4f \n\n',Q_ax,corr_ax);
fprintf(fid_out,'Quality Factor for Axially Symmetric Model = %5.3f ; Corr.coeff = %7.4f \n\n',Q_ax,corr_ax);

end


end


%=======================================================================================
%                                                                                      =
%===============================  FULLY ANISOTROPIC MODEL ==============================
%                                                                                      =
%=======================================================================================

fprintf('\n');
res=input('Fully Anisotropic model? [1]-->yes ; [0]-->no : ');
if isempty(res), res=1; end
if (res==1)
    flag_case=flag_case+1;



count=[1:1:niter];    %-- count number of synthetic data, used for drawing
start=ones(niter,6); 


%------ Array of random initial guess to search the best solution-----
for i=1:niter
    start(:,1)=guessMC_full(1)+(guessMC_full(2)-guessMC_full(1))*rand(niter,1);
    start(:,2)=guessMC_full(3)+(guessMC_full(4)-guessMC_full(3))*rand(niter,1);
    start(:,3)=guessMC_full(5)+(guessMC_full(6)-guessMC_full(5))*rand(niter,1);
    start(:,4)=guessMC_full(7)+(guessMC_full(8)-guessMC_full(7))*rand(niter,1);
    start(:,5)=guessMC_full(9)+(guessMC_full(10)-guessMC_full(9))*rand(niter,1);
    start(:,6)=guessMC_full(11)+(guessMC_full(12)-guessMC_full(11))*rand(niter,1); 
end

fprintf('Search Range\n');
fprintf('%5.1f < Dxx   < %5.1f\n',guessMC_full(1),guessMC_full(2));
fprintf('%5.1f < Dyy   < %5.1f\n',guessMC_full(3),guessMC_full(4));
fprintf('%5.1f < Dzz   < %5.1f\n',guessMC_full(5),guessMC_full(6));
fprintf('%5.1f < alpha < %5.1f\n',guessMC_full(7)*180.0/pi,guessMC_full(8)*180.0/pi);
fprintf('%5.1f < beta  < %5.1f\n',guessMC_full(9)*180.0/pi,guessMC_full(10)*180.0/pi);
fprintf('%5.1f < gamma < %5.1f\n',guessMC_full(11)*180.0/pi,guessMC_full(12)*180.0/pi);


%----------- calculate difference between calculated and exp data -----
mapendlsq=ones(niter,7);

%----- bounding conditions in the fitting procedure-------
lb=[bound_full(1) bound_full(3) bound_full(5) bound_full(7) bound_full(9) bound_full(11)];   %lower bound for Dxx Dyy Dzz alpha beta gamma
ub=[bound_full(2) bound_full(4) bound_full(6) bound_full(8) bound_full(10) bound_full(12)];   %upper bound for Dxx Dyy Dzz alpha beta gamma
  




if (alg==1),
    %options_full=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000,'LargeScale','on','LevenbergMarquardt','on');
    options_full=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000);
else
    options_full=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000);
end


fprintf('\n');    
fprintf('********************************************************************************************\n');
fprintf('****************  Starting Monte Carlo simulation with %0.f initial values   *****************\n',niter);
fprintf('****************               for fully anisotropic model                 *****************\n');
fprintf('********************************************************************************************\n');
tic

fprintf('Calculation in progress\n');
fprintf('0%%');
for i=1:niter-6
    fprintf(' ');
end
fprintf('100%% \n');


%----begin the MC search-------
for i=1:niter
     x=lsqnonlin(@rho_full2,start(i,:),lb,ub,options_full,ratio,sigma,vcoor,wNfull); %---- LVM algorithm----
     errlsq=rho_full(x,ratio,sigma,vcoor,wNfull);
     
     %------ calc progress ---------
   
    fprintf('>');
     
    
    %---- array of results with LM-----
    mapendlsq(i,1)=x(1);
    mapendlsq(i,2)=x(2);
    mapendlsq(i,3)=x(3);
    mapendlsq(i,4)=x(4);
    mapendlsq(i,5)=x(5);
    mapendlsq(i,6)=x(6);
    mapendlsq(i,7)=errlsq;
    
end

fprintf('\n');
toc



%---- plot of chi^2 ------
figure(3);
subplot(2,1,2);
hist(mapendlsq(:,7),20)
title('\fontsize{12}Fully Anisotropic model');
xlabel('\fontsize{12}\chi^{2}');
ylabel('\fontsize{12}N','Rotation',0);


%---reorder matrix of results to pick the best target function-----
%---- with the convention  Dzz>Dyy>Dxx ------------Initial valuesInitial values----------------

sortlsq_full=sortrows(mapendlsq,7);
j=1;
u=0;
lim=sortlsq_full(1,7);

for i=1:niter
    if((sortlsq_full(i,7)-lim)<=10^-5) %count number of lowest target functions
        u=u+1;
    end
end

der=ones(u,7);
for i=1:niter
    if ((sortlsq_full(i,7)-lim)<=10^-5); %build array containing only the best target functions
        der(j,:)=sortlsq_full(i,:);
        
    end
    j=j+1;
end
sol=ones(1,7);
sortder=sortrows(der,1);
for i=1:u
    if ((sortder(i,1)<sortder(i,2))&(sortder(i,2)<sortder(i,3))); %sort so as to obey convention Dzz>Dyy>Dxx
        sol(1,:)=sortder(i,:);
    
    end
end

if (((sol(1,3)-sol(1,2))<=1e-10)|(sol(1,3)-sol(1,1)<=1e-10)),
    fprintf(' \n');
    disp('warning!!!!! Unable to find a solution where Dzz>Dyy>Dxx -----> try to increase niter');
    fprintf('\n');
    fprintf(fid_out,'warning!!!!! Unable to find a solution where Dzz>Dyy>Dxx -----> try to increase niter\n\n');
    fclose(fid_out);
    return
end


%------ prepare report-------
zres=sol;

D2tau_lsq=[10^2/(6*zres(1)) 10^2/(6*zres(2)) 10^2/(6*zres(3))]; %calculates TAUx TAUy TAUz for LVM

tauclsq=(D2tc(zres(1),zres(2),zres(3)));                          %calculates TAUc

anilsq=aniso(zres(1),zres(2),zres(3));                           %Anisotropy

rhomblsq=rhomb(zres(1),zres(2),zres(3));                        %Rhombicity

[DztoDz,DztoDy,DztoDx]=rap(zres(1),zres(2),zres(3));            %calculates ratio Dzz/Dzz Dzz/Dyy Dzz/Dxx

fprintf('\n');
fprintf('=================== FULLY ANISOTROPIC MODEL ======================= \n');
fprintf('Fully Anisotropic model, Final results : \n');
fprintf('Dxx = %5.2f*10^7  ; Dyy = %5.2f*10^7 ;  Dzz = %5.2f*10^7  \n',zres(1),zres(2),zres(3));
fprintf('alpha =%4.0f deg ; beta =%4.0f deg ; gamma =%4.0f deg ; Chi2 =%8.3f ; Chi2/df=%7.3f \n',zres(4)*180.0/pi,zres(5)*180./pi,zres(6)*180./pi,zres(7),zres(7)/df_full);
fprintf('TAUx =%5.2f ns  TAUy =%5.2f ns  TAUz =%5.2f ns \n',D2tau_lsq(1),D2tau_lsq(2),D2tau_lsq(3));
fprintf('TAUc =%5.2f ns ; anisotropy =%5.2f ; rhombicity = %5.3f ; Dyy/Dxx =%6.3f  Dzz/Dxx =%6.3f \n',tauclsq,anilsq,rhomblsq,DztoDy,DztoDx);
fprintf('=================================================================== \n');




%---- save in output file -----------------
fprintf('\n\n');
fprintf(fid_out,'=================== FULLY ANISOTROPIC MODEL ===================== \n');
fprintf(fid_out,'Fully Anisotropic model, Final results : \n\n');
fprintf(fid_out,'Dxx = %5.2f*10^7   Dyy = %5.2f*10^7   Dzz = %5.2f*10^7  \n\n',zres(1),zres(2),zres(3));
fprintf(fid_out,'alpha = %4.0f deg  beta = %4.0f deg  gamma = %4.0f deg  Chi2 =%8.3f Chi2/df=%7.3f \n\n',zres(4)*180.0/pi,zres(5)*180./pi,zres(6)*180./pi,zres(7),zres(7)/df_full);
fprintf(fid_out,'TAUx =%5.2f ns  TAUy =%5.2f ns  TAUz =%5.2f ns \n\n',D2tau_lsq(1),D2tau_lsq(2),D2tau_lsq(3));
fprintf(fid_out,'TAUc =%5.2f ns ; anisotropy =%5.2f ; rhombicity = %5.3f ; Dyy/Dxx =%6.3f  Dzz/Dxx =%6.3f \n\n',tauclsq,anilsq,rhomblsq,DztoDy,DztoDx);
fprintf(fid_out,'================================================================= \n\n');



%------ draw vectors--------

rotend=(rotation_matrix(zres(4),zres(5),zres(6))).*1.5;
figure(2);
hold on
draw_vect(rotend,'g',2.5);
axis square





%------ Confidence limits ------------------------------

res=input('calculate errors ? (It can be time consuming depending on your computer speed.....)[1]-->yes ; [0]-->no : ');
if isempty(res), res=1; end
if (res==1),
    
full_par3=[zres(1),zres(2),zres(3),zres(4),zres(5),zres(6)];

matrix_full=rotation_matrix(zres(4),zres(5),zres(6));
dR_dalpha=matrix_dalpha(zres(4),zres(5),zres(6));
dR_dbeta=matrix_dbeta(zres(4),zres(5),zres(6));
dR_dgamma=matrix_dgamma(zres(4),zres(5),zres(6));

mat_full3=invcovar_fullax_angle(full_par3,sigma,matrix_full,dR_dalpha,dR_dbeta,dR_dgamma,vcoor,wNfull);
[me_full3,sd_full3]=r2r1conf(full_par3,zres(7)/df_full,mat_full3,mc_max);
%note that in me_full and sd_full the DIff tensor values are multiplied by 10^7. 

aniso2=aniso(me_full3(1),me_full3(2),me_full3(3));          %should take this out

[error_tauc,error_aniso,error_rhomb]=diff2ani_err(full_par3(1:3),sd_full3(1:3)/1e7);


fprintf('\n');
fprintf('*********** MONTE CARLO error report for fully anisotropic model ***********\n');
fprintf('Dxx = %5.2f +/- %5.2f *10^7\n',me_full3(1)/1e7,sd_full3(1)/1e7);
fprintf('Dyy = %5.2f +/- %5.2f *10^7\n',me_full3(2)/1e7,sd_full3(2)/1e7);
fprintf('Dzz = %5.2f +/- %5.2f *10^7\n',me_full3(3)/1e7,sd_full3(3)/1e7);
fprintf('****************************************************************************\n');


%df comment: In the final resport I use MC errors but the actual Dx,Dy,Dz solutions not the average ones

fprintf('\n');
fprintf('*********** FINAL REPORT for fully anisotropic model ***********\n');
fprintf('Dxx = %5.2f +/- %5.2f *10^7\n',full_par3(1),sd_full3(1)/1e7);
fprintf('Dyy = %5.2f +/- %5.2f *10^7\n',full_par3(2),sd_full3(2)/1e7);
fprintf('Dzz = %5.2f +/- %5.2f *10^7\n',full_par3(3),sd_full3(3)/1e7);
fprintf('alpha = %4.0f +/- %3.0f deg\n',full_par3(4)*180.0/pi,sd_full3(4)*180.0/pi);
fprintf('beta = %4.0f +/- %3.0f deg\n',full_par3(5)*180.0/pi,sd_full3(5)*180.0/pi);
fprintf('gamma = %4.0f +/- %3.0f deg\n',full_par3(6)*180.0/pi,sd_full3(6)*180.0/pi);
fprintf('TAUc = %5.2f +/- %5.2f ns\n',tauclsq,error_tauc);
fprintf('anisotropy = %5.2f +/- %5.2f\n',anilsq,error_aniso); 
fprintf('rhombicity = %6.3f +/- %6.3f\n',rhomblsq,error_rhomb);
fprintf('****************************************************************************\n');
fprintf('(run <euler2all.m> to find alternative equivalent Euler angles)\n\n');


%---- save in output file ------
fprintf(fid_out,'\n\n');
fprintf(fid_out,'************** MONTE CARLO error report for fully anisotropic model *************\n');
fprintf(fid_out,'Dxx = %5.2f +/- %5.2f *10^7\n',me_full3(1)/1e7,sd_full3(1)/1e7);
fprintf(fid_out,'Dyy = %5.2f +/- %5.2f *10^7\n',me_full3(2)/1e7,sd_full3(2)/1e7);
fprintf(fid_out,'Dzz = %5.2f +/- %5.2f *10^7\n',me_full3(3)/1e7,sd_full3(3)/1e7);
fprintf(fid_out,'*********************************************************************************\n\n');

fprintf(fid_out,'\n');
fprintf(fid_out,'*********** FINAL REPORT for fully anisotropic model ***********\n');
fprintf(fid_out,'Dxx = %5.2f +/- %5.2f *10^7\n',full_par3(1),sd_full3(1)/1e7);
fprintf(fid_out,'Dyy = %5.2f +/- %5.2f *10^7\n',full_par3(2),sd_full3(2)/1e7);
fprintf(fid_out,'Dzz = %5.2f +/- %5.2f *10^7\n',full_par3(3),sd_full3(3)/1e7);
fprintf(fid_out,'alpha = %4.0f +/- %3.0f deg\n',full_par3(4)*180.0/pi,sd_full3(4)*180.0/pi);
fprintf(fid_out,'beta = %4.0f +/- %4.0f deg\n',full_par3(5)*180.0/pi,sd_full3(5)*180.0/pi);
fprintf(fid_out,'gamma = %4.0f +/- %4.0f deg\n',full_par3(6)*180.0/pi,sd_full3(6)*180.0/pi);
fprintf(fid_out,'TAUc = %5.2f +/- %5.2f ns\n',tauclsq,error_tauc);
fprintf(fid_out,'anisotropy = %5.2f +/- %5.2f\n',anilsq,error_aniso); 
fprintf(fid_out,'rhombicity = %6.3f +/- %6.3f\n',rhomblsq,error_rhomb);
fprintf(fid_out,'****************************************************************************\n\n');
fprintf(fid_out,'(run <euler2all.m> to find alternative equivalent Euler angles)\n');


end
%----- statistical significance using F-test ------------------
if (flag_case==2)

%----- full and axially symmetric -----------------------------

disp(' ');
[P_fa,F_fa]=ftestcmp(zres(7)/df_full,df_full,zres_ax(5)/df,df);
disp(['F-STATISTICS FULL. ANISO VS. AXIALLY-SYMMETRIC : F= ',num2str(F_fa),'  P= ',num2str(P_fa)]);



%---- save in output file --------------------
fprintf(fid_out,'\n\n');
fprintf(fid_out,'F-STATISTICS FULL. ANISO VS. AXIALLY-SYMMETRIC : F= %5.2g ; P= %5.2E\n',F_fa,P_fa);


%---- full and isotropic --------------------------------------
disp(' ');
[P_fi,F_fi]=ftestcmp(zres(7)/df_full,df_full,chi2iso/dfiso,dfiso);
disp(['F-STATISTICS FULL. ANISO VS. ISOTROPIC : F= ',num2str(F_fi),'  P= ',num2str(P_fi)]);
fprintf('\n');

%---- save in output file --------------------
fprintf(fid_out,'\n\n');
fprintf(fid_out,'F-STATISTICS FULL. ANISO VS. ISOTROPIC : F= %5.2g ; P= %5.2E\n\n',F_fi,P_fi);

end


%----- visualise points ---------------------

res=input('visualize points? [1]-->yes ; [0]-->no : ');
if isempty(res), res=1; end

if (res==1),
    if ((flag_case==2)&(flag_visu==1)),
[Q_full,fan_calc,corr_full]=plot_orientation(zres(1),zres(2),zres(3),zres(4),zres(5),zres(6),freq,r2,r1,NH,reslist,r3,wNfull,l_lim,u_lim,flag_case);
%end

%if ((flag_case==2)&(flag_visu~=1)),
%plot_orientation(zres(1),zres(2),zres(3),zres(4),zres(5),zres(6),freq,r2,r1,NH,reslist,r3,wNfull,0,0,1);
    else
[Q_full,fan_calc,corr_full]=plot_orientation(zres(1),zres(2),zres(3),zres(4),zres(5),zres(6),freq,r2,r1,NH,reslist,r3,wNfull,0,0,0); 

    end
fprintf('\n');
fprintf('Quality Factor for Fully Anisotropic Model = %5.3f ; Corr.coeff = %7.4f \n',Q_full,corr_full);
fprintf(fid_out,'Quality Factor for Fully Anisotropic Model = %5.3f ; Corr.coeff = %7.4f \n\n',Q_full,corr_full);
end
end

tab_gen=[init_tab,fan_calc,ax_calc];
mat2ascii('tab_gen.txt',tab_gen);   %record in file
fprintf('\n');
fclose(fid_out);
return

%===================================================================






