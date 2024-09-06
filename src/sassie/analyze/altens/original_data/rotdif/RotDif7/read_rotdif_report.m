function z=read_rotdif_report(fname,dirname);
%-------------------------------------------
%   df-mar-04
%   read an ascii file -- output of Rotdif
%   and extract fit results
%-------------------------------------------
if nargin < 2, dirname='./'; end    %default
z=[];
%dirname='./VEAN/CORE/';
fullfname=fullfile(dirname,fname);
txt=readasci(fullfname);
nlin=size(txt,1);
for ii=1:nlin,
    if (~isempty(findstr('alpha',txt(ii,:)))  & ~isempty(findstr('beta',txt(ii,:))) & ~isempty(findstr('gamma',txt(ii,:)))),
        tmp = txt(ii,:);
        ind_eq = findstr('=',tmp);
        ind_deg = findstr('deg',tmp);
        ind_Chi = findstr('Chi',tmp);
        alpha = str2num(tmp(ind_eq(1)+1:ind_deg(1)-1));    
        beta = str2num(tmp(ind_eq(2)+1:ind_deg(2)-1));   
        gamma = str2num(tmp(ind_eq(3)+1:ind_deg(3)-1));
        Chi2_df = str2num(tmp(ind_eq(4)+1:ind_Chi(2)-1));
        Chi2 = str2num(tmp(ind_eq(5)+1:end));
        z=[alpha,beta,gamma,Chi2_df,Chi2];
        return
    end
end
