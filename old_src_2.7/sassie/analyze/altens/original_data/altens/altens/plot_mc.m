function plot_mc(ang,param)

figure(4)
clf
%[n1,xout1]=hist(ang(:,1));
subplot(2,3,1)
[N,Xout]=hist(ang(:,1),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{16}\alpha');


%[n2,xout2]=hist(ang(:,2));
hold on
subplot(2,3,2)
[N,Xout]=hist(ang(:,2),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{16}\beta');



%[n3,xout3]=hist(ang(:,3));
hold on
subplot(2,3,3)
[N,Xout]=hist(ang(:,3),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{16}\gamma');


%[n4,xout4]=hist(param(:,1));
hold on
subplot(2,3,4)
[N,Xout]=hist(param(:,1),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{14} Sxx');


%[n5,xout5]=hist(param(:,2));
hold on
subplot(2,3,5)
[N,Xout]=hist(param(:,2),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{14} Syy');


%[n6,xout6]=hist(param(:,3));
hold on
subplot(2,3,6)
[N,Xout]=hist(param(:,3),20);
h=bar(Xout,N);
set(h,'FaceColor','r');
xlabel('\fontsize{14} Szz');
