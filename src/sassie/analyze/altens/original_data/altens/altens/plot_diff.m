function plot_diff(new,old,R)
if nargin < 3, R=[]; end
[old,new];
clf
grid on
plot(old,new,'.r','MarkerSize',20);




upper_bound=max(old);
lower_bound=min(old);

hold on
co1=corrcoef(old,new);
xlabel('\fontsize{14}experimental RDC');
ylabel('\fontsize{14}calculated RDC');
axis([lower_bound*2 upper_bound*2 lower_bound*2 upper_bound*2]);

plot(linspace(lower_bound*2,upper_bound*2,40),linspace(lower_bound*2,upper_bound*2,40));

title(['\fontsize{14} corr. coeff. = ',num2str(co1(1,2)),'  Rfactor =',num2str(R)]);