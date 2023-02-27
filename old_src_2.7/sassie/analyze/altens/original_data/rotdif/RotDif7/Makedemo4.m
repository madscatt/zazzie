%makedemo
freq=600.13;
TAUc=8.0;
Dx=1.7;
Dy=2.1;
Dz=2.3;



tab=[];
tab_param=[];
tab_final_param=[];
        


for i=1:20
    i
    tab_Dy=ones(50,1).*(Dy/Dx)
    [thetamap,ratio_th1,ratio_th2]=plot_orientation_demo(Dx,Dy,Dz);
    tab_final_param=[tab_final_param,[thetamap,tab_Dy,ratio_th1,ratio_th2]]
    pause
    Dy=Dy-0.02;
    
end
tab_final_param
