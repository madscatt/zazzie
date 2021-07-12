function [X,Y,matrix]=createplot2(tab,col)


x=0.005:0.005:0.1
y=tab(1:10,4)'

matrix=reshape(tab(:,col),10,20)

figure
surf(x,y,matrix)
shading interp
%alpha 0.9
axis square



