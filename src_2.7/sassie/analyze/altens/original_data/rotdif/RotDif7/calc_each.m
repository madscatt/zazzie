function t_each=calc_each(omeg,r)


a=(1./r).*(3/4);
b=sqrt(a-1);
c=(1/omeg).*b;

t_each=c;
return