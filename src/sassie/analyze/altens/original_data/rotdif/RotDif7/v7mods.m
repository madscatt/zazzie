two changes to lsqncommon.m
ln 198 rep %Jstr = optimget(options,'JacobPattern',[]);
        Jstr=[];
        
        ln 52 %JH change 8/06 for Matlab v7 %%%%%%%%%%%%%%%%
%mtxmpy = optimget(options,'JacobMult',[]); % use old
mtxmpy=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snls
 ln 90 rep %active_tol = optimget(options,'ActiveConstrTol',sqrt(eps)); % leave old optimget
active_tol=sqrt(eps);
ln 95 %pcmtx = optimget(options,'Preconditioner','aprecon') ; % leave old optimget
pcmtx ='aprecon';

ln 139 
%showstat = optimget(options,'showstatus','off'); % no default, so leave with slow optimget
showstat='off';

changed sfdnls to xsfdnls 
change function call
ln 183 snls.m
ln 331 snls.m


r2r1visu.m ln 107
rep %bar(reslist,diff,'ro');
bar('v6',reslist,diff,'ro');

plot_orientation at 141 rep %bar(rlist,diff,'ro')
bar('v6',rlist,diff,'ro')

rotdif 261
%options=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000,'LargeScale','on','LevenbergMarquardt','on');
    options=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000);
    625
%options_full=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000,'LargeScale','on','LevenbergMarquardt','on');
    options_full=optimset('TolX',1.e-15,'TolFun',1e-15,'MaxIter',8000,'maxFunEvals',20000);

    
