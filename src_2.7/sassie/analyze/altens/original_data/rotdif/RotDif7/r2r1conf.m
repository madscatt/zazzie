function [me_mc,sd_mc]=r2r1conf(PAR0,VARmin,invCov,mc_max)
%-----------------------------------------------------------------------------
% df-dec-98	
%       estimate errors in the microdynamic params
%       by MC simulation using the inverse covariance
%	matrix approach and chi2-boundaries
% 	INPUT:  PAR0 - row-vector of optim. params
%          VARmin - chi2(@par0)/df variation at min
%          invCov - inverse Covar. mtrx
%          mc_max - max number of MC points
%          mcmode=2 - toggles the visualization mode
% Copyright (c) 1996-97 by David Fushman, The Rockefeller University
%-----------------------------------------------------------------------------
 dchi2=[1 2.3 3.53 4.72 5.89 7.04];   % 68.3% conf.level
 %dchi2=[2.71 4.61 6.25 7.78 9.24 10.6];
 npar=length(PAR0);
 PAR0(1)=PAR0(1)*1e7;
 PAR0(2)=PAR0(2)*1e7;
 PAR0(3)=PAR0(3)*1e7;
 me_mc=zeros(1,npar);
 sd_mc=zeros(1,npar);
 chi2bound=dchi2(npar)*VARmin;
 if VARmin<0.01,
   d0=0.0002*ones(1,npar);                        	%starting deviations
 else
   d0=0.02*ones(1,npar);
 end

 sim=zeros(mc_max,npar);
 mc_i=1;
                                     %visualization = off
   RES_mc=zeros(2,npar);
   REStmp=zeros(3,npar);
   while mc_i<=mc_max,                          %MC-simulation
     sim=2*(rand(mc_max,npar)-0.5);  %simulate data in the [-1,1] interval
    
     for itmp=1:mc_max,
        dd=d0.*sim(itmp,:);
        dPAR=PAR0.*dd;
      
        if dPAR*invCov*dPAR' <= chi2bound;
          border=ones(1,npar);
         
          for iii=1:npar,
            if abs(dd(iii))>0.95*d0(iii),       %crossing the border!!!
                border(iii)=2;
                
            end
          end
          if border(:)==1,                      %border control:ok
             REStmp(1:2,:)=RES_mc;
             REStmp(3,:)=dPAR;
             RES_mc(1,:)=max(REStmp);
             RES_mc(2,:)=min(REStmp);
             mc_i=mc_i+1;
         
          else                                  %reached the border
            d0=d0.*border;
            if d0(:)<=5.121,	                %or 0.02 & 2.561(240%)
              mc_i=1;                           %back to the beginning
            else
              for iipar=1:npar,
                if d0(iipar)>5.121,
                  disp(['!!! par(',num2str(iipar),') error > 500% !!!']);
                  RES_mc(:,iipar)=NaN*ones(2,1);   %mark as NaN's
                end
              end
              mc_i=mc_max+1;
              break
            end
          end
        end
     end                                        %for itmp=1:mc_max,
   end                                          %MC-simulation
   me_mc=(RES_mc(1,:)+RES_mc(2,:))/2;
   sd_mc=(RES_mc(1,:)-RES_mc(2,:))/2;
  me_mc=me_mc+PAR0;				%adjustment
return
%============================================================

