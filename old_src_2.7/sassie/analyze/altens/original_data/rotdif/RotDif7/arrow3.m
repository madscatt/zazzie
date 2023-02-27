function H=arrow3(p1,p2,s,w,h,ip)
% ARROW3
%   ARROW3(P1,P2) will draw vector lines (2D/3D) from P1 to P2 with
%   arrowheads, where P1 and P2 can be either nx2 matrices (for 2D),
%   or nx3 matrices (for 3D).
%
%   ARROW3(P1,P2,S,W,H,IP) can be used to specify properties of the
%   line and arrowhead.  S is a character string made with one element
%   from any or all of the following 3 columns:
%
%      Color Switches        LineStyle             LineWidth
%      -------------------   -------------------   --------------------
%      k   black (default)   -   solid (default)   0.5 points (default)
%      y   yellow            :   dotted            0   no lines
%      m   magenta           -.  dashdot           /   LineWidthOrder
%      c   cyan              --  dashed
%      r   red               *   LineStyleOrder
%      g   green
%      b   blue
%      w   white
%      x   random named color
%      o   ColorOrder
%
%   LineWidth units are points (1/72 inch); use 0 to omit lines.
%   Invalid characters in S will be ignored and replaced by default
%   settings.  W is a vector of arrowhead widths.  H (default = 3W) is
%   a vector of arrowhead heights.  If a value is provided for IP, an
%   initial point marker will be plotted with width IP; use IP = 0 for
%   default width W.  For linear plots with equal axes, the units of W,
%   H, and IP are the same as those of the coordinate data (P1,P2).
%
%   Note that MATLAB cycles through the line styles defined by the
%   LineStyleOrder property only after using all colors defined by
%   the ColorOrder property.  If however, the global variable
%   LineWidthOrder is defined, and LineWidth is specified with '/',
%   then each line will be drawn with sequential color, linestyle,
%   and linewidth.
%
%   Plotting lines with a single color, linestyle, and linewidth is
%   faster than plotting lines with multiple colors and/or linestyles.
%   Plotting lines with multiple linewidths is slower still.
%
%   If a particular aspect ratio is required, use DASPECT or PBASPECT
%   before calling ARROW3.  Changing DASPECT or PBASPECT after calling
%   ARROW3 may change the appearance of arrowheads and initial point
%   markers.  Arrow3 sets DataAspectRatioMode to manual for linear
%   plots and sets PlotBoxAspectRatioMode to manual for log plots.
%
%   Usage Examples:
%
%       % 2D vectors
%       Arrow3([0 0],[1 3])
%       Arrow3([0 0],[1 2],'-.b')
%       Arrow3([0 0],[10 10],'--x2',1)
%       Arrow3(zeros(10,2),50*rand(10,2),'x',1,3)
%       Arrow3(zeros(10,2),[10*rand(10,1),500*rand(10,1)],'r')
%       Arrow3(10*rand(10,2),50*rand(10,2),'x',1,[],1)
%
%       % 3D vectors
%       Arrow3([0 0 0],[1 1 1])
%       Arrow3(zeros(20,3),50*rand(20,3),'--x1.5',2)
%       Arrow3(zeros(100,3),50*rand(100,3),'x',1,3)
%       Arrow3(zeros(10,3),[10*rand(10,1),500*rand(10,1),50*rand(10,1)],'b')
%       Arrow3(10*rand(10,3),50*rand(10,3),'x',[],[],0)
%
%       % Just for fun
%       Arrow3(zeros(100,3),50*rand(100,3),'x',10,3);
%       light('Position',[-10 -10 -10],'Style','local');
%       light('Position',[60,60,60]), lighting gouraud; alpha(.95)
%
%       % ColorOrder, LineStyleOrder, and LineWidthOrder
%       set(gca,'LineStyleOrder',{'-','--','-.',':'})
%       set(gca,'ColorOrder',[1,0,0;0,0,1;0.25,0.75,0.25;0,0,0])
%       global LineWidthOrder, LineWidthOrder=[1,2,4,8];
%       w=[5,10,15,20]; h=[20,40,30,40];
%       Arrow3(zeros(4,2),[10*rand(4,1),500*rand(4,1)],'o*/',w,h,10)
%
%       % Log plot
%       loglog([1e2,1e8],[1e-2,1e-1],'wo','MarkerSize',eps), grid on, hold on
%       p1=repmat([1e3,2e-2],15,1);
%       q1=[1e7,1e6,1e5,1e4,1e3,1e7,1e7,1e7,1e7,1e7,1e7,1e6,1e5,1e4,1e3];
%       q2=1e-2*[9,9,9,9,9,7,5,4,3,2,1,1,1,1,1]; p2=[q1',q2'];
%       set(gca,'ColorOrder',rand(15,3));
%       Arrow3(p1,p2,'o'), hold off
%
%   Copyright(c)2002, Version 4.00
%       Jeff Chang <cpmame@hotmail.com>
%       Tom Davis  <tdavis@eng.usf.edu>

%   Revision History:
%
%       10/12/02   - Added global LineWidthOrder,
%                    vectorized W, H and IP (TD)
%       10/05/02   - Changed CLF to CLA for subplot support,
%                    added ColorOrder and LineStyleOrder support (TD)
%       04/27/02   - Minor log plot revisions (TD)
%       03/26/02   - Added log plot support (TD)
%       03/24/02   - Adaptive grid spacing control to trade off
%                    appearance vs. speed based on size of matrix (JC)
%       03/16/02   - Added "axis tight" for improved appearance (JC)
%       03/12/02   - Added initial point marker (TD)
%       03/03/02   - Added aspect ratio support (TD)
%       03/02/02   - Enchance program's user friendliness (JC)
%                    (lump Color, LineStyle, and LineWidth together)
%       03/01/02   - Replaced call to ROTATE (TD)
%       02/28/02   - Modified line plotting,
%                    added linewidth and linestyle (TD)
%       02/27/02   - Minor enhancements on 3D appearance (JC)
%       02/26/02   - Minor enhancements for speed (TD&JC)
%       02/26/02   - Optimise PLOT3 and SURF for speed (TD)
%       02/25/02   - Return handler, error handling, color effect,
%                    generalize for 2D/3D vectors (JC)
%       02/24/02   - Optimise PLOT3 and SURF for speed (TD)
%       02/23/02   - First release (JC&TD)
%

%===============================================================
% Error Checking
if ishold, restore=1; else restore=0; end
if nargin<2
    error([upper(mfilename),' requires at least two input arguments'])
end
[r1,c1]=size(p1); [r2,c2]=size(p2); n=r1;
if r1~=r2, error('P1 and P2 must have same number of rows'), end
if c1~=c2, error('P1 and P2 must have same number of columns'), end
if c1==1
    error(['Invalid input, type HELP ',upper(mfilename),' for usage examples'])
end
if c1==2, p1=[p1,zeros(n,1)]; p2=[p2,zeros(n,1)]; end
if ~restore, cla, hold on, xys=0; end, view(c1), F=gca;
L=get(F,'LineStyleOrder'); C=get(F,'ColorOrder');
if restore
    xs=strcmp(get(F,'xscale'),'log');
    ys=strcmp(get(F,'yscale'),'log');
    zs=strcmp(get(F,'zscale'),'log');
    if zs, error('Z log scale not supported'), end
    xys=xs+ys;
    if xys & c1==3, error('3D log plot not supported'), end
end

%===============================================================
% Style Control
vc=['ymcrgbkwxo'];                          % valid color codes
cn=[1,1,0;1,0,1;0,1,1;1,0,0;0,1,0;0,0,1;0,0,0;1,1,1];
if nargin<3, [c,ls,lw]=LocalValidateCLSW;   % default Color, LineStyle/Width
else, 
    [c,ls,lw]=LocalValidateCLSW(s);
    if c=='x'                               % random named color (less white)
        c=cn(rem(randperm(n),7)+1,:);
        set(F,'ColorOrder',c)
    elseif c=='o'                           % use default ColorOrder
        c=repmat(C,ceil(n/size(C,1)),1);
    elseif ~sum(vc==c),
        warning('Invalid color switch, default color (black) will be used');
        c='k';
    end
end
if ls~='*', set(F,'LineStyleOrder',ls), end
if length(c)==1, set(F,'ColorOrder',cn(find(vc==c),:)), end
if lw=='/', global LineWidthOrder,
    if length(LineWidthOrder)
        lw=repmat(LineWidthOrder(:),ceil(n/length(LineWidthOrder)),1);
    else
        warning('Undefined LineWidthOrder, default width (0.5) will be used');
        lw=0.5;
    end
end

%===============================================================
% Log Plot
if xys
    axis auto
    set(F,'units','points')
    pos=get(F,'position');
    set(F,'units','normalized')
    xr=get(F,'xlim'); yr=get(F,'ylim');
    if strcmp(get(F,'PlotBoxAspectRatioMode'),'auto')
        set(F,'PlotBoxAspectRatio',[pos(3),pos(4),1])
    end
    par=get(F,'PlotBoxAspectRatio');
    set(F,'DataAspectRatio',[par(2),par(1),par(3)])
    % map coordinates onto unit square
    q=[p1;p2];
    if xs, xr=log10(xr); q(:,1)=log10(q(:,1)); end
    if ys, yr=log10(yr); q(:,2)=log10(q(:,2)); end
    q=q-repmat([xr(1),yr(1),0],2*n,1);
    dx=xr(2)-xr(1); dy=yr(2)-yr(1);
    q=q*diag([1/dx,1/dy,1]);
    q1=q(1:n,:); q2=q(n+1:end,:);
end

%===============================================================
% Line
if length(lw)==1
    if lw>0
        if length(c)==1 & ls~='*' % single color, linestyle, and linewidth
            P=zeros(3*n,3); i=1:n;
            P(3*i-2,:)=p1(i,:); P(3*i-1,:)=p2(i,:); P(3*i,1)=NaN;
            H1=plot3(P(:,1),P(:,2),P(:,3),'LineWidth',lw,'Tag','arrow3');
        else                      % single linewidth
            H1=plot3([p1(:,1),p2(:,1)]',[p1(:,2),p2(:,2)]',...
                [p1(:,3),p2(:,3)]','LineWidth',lw,'Tag','arrow3');
        end
    else, H1=[];
    end
    if length(c)==1, c=repmat(c,n,1); end
else                              % use LineWidthOrder
    if length(c)==1, c=repmat(c,n,1); end
    ls=repmat(cellstr(L),ceil(n/size(L,1)),1);
    H1=zeros(n,1);
    for i=1:n
        H1(i)=plot3([p1(i,1),p2(i,1)],[p1(i,2),p2(i,2)],[p1(i,3),p2(i,3)], ...
            ls{i},'Color',c(i,:),'LineWidth',lw(i),'Tag','arrow3');
    end
end

%===============================================================
% Scale
ar=get(F,'DataAspectRatio'); ar=sqrt(3)*ar/norm(ar);
set(F,'DataAspectRatioMode','manual')
if nargin<4 | isempty(w)              % width
    if xys, w=1;
    else
        xr=get(F,'xlim'); yr=get(F,'ylim'); zr=get(F,'zlim');
        w=norm([xr(2)-xr(1),yr(2)-yr(1),zr(2)-zr(1)])/72;
    end
end
if xys, w=w/72; end
if nargin<5 | isempty(h), h=3*w; end  % height
w=repmat(w(:),ceil(n/length(w)),1);
h=repmat(h(:),ceil(n/length(h)),1);

%===============================================================
% Arrowhead
if xys, p1=q1; p2=q2; end
W=(p1-p2)./repmat(ar,n,1);            % new z direction
W=W./repmat(sqrt(sum(W.*W,2)),1,3);   % unit vector
U=[-W(:,2),W(:,1),zeros(n,1)];        % new x direction
N=sqrt(sum(U.*U,2));                  % norm(U)
i=find(N<eps); j=length(i);
U(i,:)=repmat([1,0,0],j,1); N(i)=ones(j,1);
U=U./repmat(N,1,3);                   % unit vector
V=cross(W,U,2);                       % new y direction

m1=30; m2=10; num=200;                % max, min grid spacing, and threshold
if n<num, m=round((m2-m1)*n/num+m1);  % adjust grid spacing automatically
else m=m2; end                        % for speed when matrix size > num

[x,y,z]=cylinder([0,1],m);
G=surf(x/2,y/2,z);
X=get(G,'XData'); Y=get(G,'YData'); Z=get(G,'ZData');
[j,k]=size(X);
H2=zeros(n,1);
for i=1:n                             % translate, rotate, and scale
    H2(i)=copyobj(G,F);
    newxyz=[w(i)*X(:),w(i)*Y(:),h(i)*Z(:)]*[U(i,:);V(i,:);W(i,:)]*...
        diag(ar)+repmat(p2(i,:),j*k,1);
    newx=reshape(newxyz(:,1),j,k);
    newy=reshape(newxyz(:,2),j,k);
    newz=reshape(newxyz(:,3),j,k);
    if xys
        newx=newx*dx+xr(1); newy=newy*dy+yr(1);
        if xs, newx=10.^newx; end
        if ys, newy=10.^newy; end
    end
    set(H2(i),'XData',newx,'YData',newy,'ZData',newz,...
        'FaceColor',c(i,:),'EdgeColor','None','Tag','arrow3')
end
delete(G)

%===============================================================
% Initial Point Marker
if nargin>5 & ~isempty(ip)
    ip=repmat(ip(:),ceil(n/length(ip)),1); if xys, ip=ip/72; end
    [x,y,z]=sphere(m);
    G=surf(x*ar(1)/2,y*ar(2)/2,z*ar(3)/2);
    X=get(G,'XData'); Y=get(G,'YData'); Z=get(G,'ZData');
    H3=zeros(n,1);
    for i=1:n                         % translate
        H3(i)=copyobj(G,F);
        if ip(i)<=0, ip(i)=w(i); end
        newx=p1(i,1)+X*ip(i); newy=p1(i,2)+Y*ip(i); newz=p1(i,3)+Z*ip(i);
        if xys
            newx=newx*dx+xr(1); newy=newy*dy+yr(1);
            if xs, newx=10.^newx; end
            if ys, newy=10.^newy; end
        end
        set(H3(i),'XData',newx,'YData',newy,'ZData',newz,...
            'FaceColor',c(i,:),'EdgeColor','None','Tag','arrow3')
    end
    delete(G)
else, H3=[];
end

%===============================================================
% Finish
if ~restore, hold off, end
set(F,'LineStyleOrder',L), set(F,'ColorOrder',C)
if c1==3, axis vis3d, end             % good for 3D view
if xys, set(F,'DataAspectRatioMode','auto')
else, set(gcf,'Renderer','OpenGL'), axis tight
end
if nargout, H=[H1(:);H2(:);H3(:)]; end

%===============================================================
% Generate valid value for color, linestyle and linewidth
function [c,ls,lw]=LocalValidateCLSW(s)
if nargin<1, c='k'; ls='-'; lw=0.5;
else
    % identify linestyle
    if findstr(s,'--'), ls='--'; s=strrep(s,'--','');
    elseif findstr(s,'-.'), ls='-.'; s=strrep(s,'-.','');
    elseif findstr(s,'-'), ls='-'; s=strrep(s,'-','');
    elseif findstr(s,':'), ls=':'; s=strrep(s,':','');
    elseif findstr(s,'*'), ls='*'; s=strrep(s,'*','');
    else, ls='-';
    end

    % identify linewidth
    tmp=double(s);
    tmp=find(tmp>45 & tmp<57);
    if length(tmp)
        if any(s(tmp)=='/'), lw='/'; else, lw=str2num(s(tmp)); end
        s(tmp)='';
    else, lw=0.5;
    end

    % identify color
    if length(s), c=s(1); else, c='k'; end
end