function [err_tau,err_ani,err_rhomb]=diff2ani_err(Dxyz,errDxyz)
%------------------------------------------
%   df-sep-15
%   convert diff tensor into tau, ani, rhomb and calculate
%   errors using error propagation equations
%   (my notebook, pp.63-64)
%------------------------------------------
      Dx=Dxyz(1);
      Dy=Dxyz(2);
      Dz=Dxyz(3);
      tau=1/(2*sum([Dx,Dy,Dz]))*1e2;
      err_tau=tau/sum([Dx,Dy,Dz])*sqrt(errDxyz*errDxyz');
      if abs(Dx-Dy) <= abs(Dy-Dz),  % prolate case
          ani=2*Dz/(Dx+Dy); 
          der_ani=[1/(Dx+Dy);1/(Dx+Dy);1/Dz];
          err_ani= ani*sqrt((errDxyz.^2)*(der_ani.^2));
          
          rhomb=3*(Dy-Dx)/(2*Dz-Dx-Dy);
          if rhomb~=0, 
            der_rhomb=[(Dz-Dy)/(abs(Dy-Dx)+eps);(Dz-Dx)/(abs(Dy-Dx)+eps);1];
            err_rhomb= 2*rhomb/(abs(2*Dz-Dx-Dy)+eps)*sqrt((errDxyz.^2)*(der_rhomb.^2));
          else
            der_rhomb=[(Dz-Dy);(Dz-Dx);0];
            err_rhomb= 2*3/(abs(2*Dz-Dx-Dy)+eps)^2*sqrt((errDxyz.^2)*(der_rhomb.^2));              
          end
          z=[tau err_tau ani rhomb err_ani err_rhomb];
      else                          % oblate case
          ani=2*Dx/(Dy+Dz);
          der_ani=[1/Dx; 1/(Dz+Dy); 1/(Dz+Dy)];
          err_ani= ani*sqrt((errDxyz.^2)*(der_ani.^2));
          
          rhomb=3*(Dy-Dz)/(2*Dx-Dy-Dz);
          if rhomb~=0, 
            der_rhomb=[1;(Dz-Dx)/(abs(Dz-Dy)+eps);(Dy-Dx)/(abs(Dz-Dy)+eps)];
            err_rhomb= 2*rhomb/(abs(2*Dx-Dy-Dz)+eps)*sqrt((errDxyz.^2)*(der_rhomb.^2));
          else
            der_rhomb=[0;(Dz-Dx);(Dy-Dx)];
            err_rhomb= 2*3/(abs(2*Dx-Dy-Dz)+eps)^2*sqrt((errDxyz.^2)*(der_rhomb.^2));              
          end
          z=[tau err_tau ani rhomb err_ani err_rhomb];
      end  

%===========================================