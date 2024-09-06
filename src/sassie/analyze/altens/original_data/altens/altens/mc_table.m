function z=mc_table(table)

num=length(table(:,1));

%------Generate random errors for RDCs---------

randn('state',sum(100*clock));
D_err=randn(num,1);


%------Create synthetic dataset incorporating errors---

z(:,1)=table(:,1)+D_err*2;
%z(:,2)=D_err.*table(:,2)+table(:,1);


return

%======================================================