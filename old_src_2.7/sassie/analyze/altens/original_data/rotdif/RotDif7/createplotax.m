function [X,Y,matrix]=createplotax(tab,col)


x=0.01:0.005:0.105
y=tab(1:10,1)'

matrix=reshape(tab(:,col),10,20)

figure
surf(x,y,matrix)
shading interp
alpha 0.9
axis square