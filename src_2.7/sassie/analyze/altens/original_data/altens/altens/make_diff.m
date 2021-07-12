function [diff]=make_diff(mat1,mat2);
% ow-july 2002 ; University of Maryland
% calculates difference between two sets of dipolar coupling
% each set corresponding to one of the doublet (J+d;J-d)
% remove peaks which are not defined



l=size(mat1,1);
diff=[];
for ii=1:l
    if ((mat1(ii,4)~=1)&(mat2(ii,4)~=1))
        diff=[diff;[mat1(ii,1) (mat1(ii,2)-mat2(ii,2))*60.81]];
    end
end

return