function mat2ascii(filename, X, fmt_string)
%MAT2ASCII - Outputs matrix X to an ASCII file
%
%Syntax:  mat2ascii(filename, X, fmt_string);
%
%Inputs:  FILENAME - Output filename
%         X - a 2-D double array
%         FMT_STRING - format string
%
%Output:  none
%
%Examples:  
%           mat2ascii('b_est.dat', b_est, '  %11.8f\t');
%           mat2ascii('stats_est.dat', stats_est, '  %11.8f\t')
%
%           X = rand(20000,10);
%           mat2ascii('foo1.dat', X, '%4.2f ')
%


if nargin < 3 , fmt_string = '%f ';  end
if nargin < 2
    error('Two input arguments required:  filename and X')
end

[nrows, ncols]  = size(X);

kpercent = strfind(fmt_string,'%');

fid = fopen(filename, 'wt');
if kpercent == 1
   fprintf(fid,[repmat(fmt_string,[1,ncols]) '\n'], X');
else
   fprintf(fid,[fmt_string '\n'], X');
end

fclose(fid);