function [res_dip_coup]=res_dip(aniso_sup,aniso_inf,iso_sup,iso_inf)
%-------------------------------------------------------
%  ow-july 2002 ; University of Maryland (modified/corrected version)
%  prepare working set of data for residual dipolar couplings
%  peaks which are not defined are eliminated
%  to be used for tensor alignment.
%  aniso_sup and aniso_inf must have the eame size
%  iso_sup and iso_inf must have the same size
%------------------------------------------------------- 




max_res=max(aniso_sup(:,1));

%------------ J+D for anisotropic phase ---------------
diff_aniso=make_diff(aniso_sup,aniso_inf);


%----------- J for isotropic phqse --------------------
diff_iso=make_diff(iso_sup,iso_inf);



%---------- construct a line at 94 Hz -----------------
line=0:max_res;
line=line'
line2=ones(max_res+1,1).*94
line=[line,line2]
res=0:max_res;



%----------- plot J+D(anisotrop) and J(isotrop) ---------
%----------- eliminate peaks which are not defined ------

max_chem=max(diff_aniso(:,2));
min_chem=min(diff_aniso(:,2));


figure(1)
clf
subplot(2,1,1)
stem(diff_aniso(:,1),diff_aniso(:,2),'fill','r')
grid on
axis([0 max_res min_chem-10 max_chem+10])
title('\fontsize{12}Anisotropic phase (J+D)');
xlabel('\fontsize{12}Residue number');
ylabel('\fontsize{12}J+D (Hz)');
hold on
plot(line(:,1),line(:,2));



subplot(2,1,2)
stem(diff_iso(:,1),diff_iso(:,2),'fill','b')
grid on
axis([0 max_res min_chem-10 max_chem+10])
title('\fontsize{12}Isotropic phase (J)');
xlabel('\fontsize{12}Residue number');
ylabel('\fontsize{12}J (Hz)');
hold on
plot(line(:,1),line(:,2),'r');





%---- calculate D ;remove peacks which are not defined for each files ----------

l=size(aniso_sup,1);
new_peak=[];
for ii=1:l
    if ((aniso_sup(ii,4)~=1)&(aniso_inf(ii,4)~=1)&(iso_sup(ii,4)~=1)&(iso_inf(ii,4)~=1))
        new_peak=[new_peak;[aniso_sup(ii,1),(((aniso_sup(ii,2)-aniso_inf(ii,2))*60.81)-((iso_sup(ii,2)-iso_inf(ii,2))*60.81))]];
    end
end

res_dip_coup=new_peak;

%---- plot D -------------------------

figure(2)
stem(res_dip_coup(:,1),res_dip_coup(:,2),'fill','g')
grid on
title('\fontsize{12}Residual Dipolar Couplings (Hz)');
xlabel('\fontsize{12}Residue number');
ylabel('\fontsize{12}D (Hz)');



return



