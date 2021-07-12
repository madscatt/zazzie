function [X,Y,matrix]=createplot(tab,col)

for i=1:10
matrix(1,i)=tab(i,col);
end

for i=2:10
    k=((i-1)*10)+1;
    for j=1:10
        matrix(i,j)=tab(k,col);
        k=k+1;
    end
end

X=0.01:0.01:0.1;
X=X';
Y=1.2:-0.02:1.02;
Y=Y';

matrix