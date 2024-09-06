function tabellip=fitpeack(rr2,dirname,xpks0)

npks=size(xpks0,1);
pick=xpks0;
%---- get parameters ---------------------------------
file2rr=[dirname,'/2rr'];
fileprocs=[dirname,'/procs'];
fileproc2s=[dirname,'/proc2s'];
%-------------------read */proc?s---------------------
[NC_proc2,offs2,sf2,si2,sw2,xdim2]=getpar(fileprocs);
[NC_proc1,offs1,sf1,si1,sw1,xdim1]=getpar(fileproc2s);



%------- spectrum limits -----------------

f1=offs1-(sw1/sf1/(si1-1))*[0:si1-1];	%f1-axis
f2=offs2-(sw2/sf2/(si2-1))*[0:si2-1];	%f2-axis

%---------------reverse axis to be consistent with spc plot--------------

f1=f1(1,end:-1:1);
f2=f2(1,end:-1:1);

%------- divide levels -------
data=rr2./100000;

%------- get max of for contour plot-------
max_cont=max(max(data));
clf
contour(f2,f1,data,round(linspace(0,100,20)));

%------ fit ellipsoid for each peack --------
tab_ellip=[];
for ipk=1:npks,
    
    
    
    range_x=[pick(ipk,3)-0.25 pick(ipk,3)+0.25]
    range_y=[pick(ipk,2)-0.1 pick(ipk,2)+0.1]
    
    [m_x,m_y]=cutlevels2(f2,f1,data,10,range_x,range_y);
    disp(['peack ',num2str(ipk),' ',num2str(m_y),' ',num2str(m_x)]);
    tab_ellip=[tab_ellip;[pick(ipk,1),m_y,m_x]];
end

tab_ellip
return
    







