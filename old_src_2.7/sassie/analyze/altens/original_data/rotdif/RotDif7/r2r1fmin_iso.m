function [xf,fvm] = r2r1fmin_iso(ax,bx,options,ratio,sigma,wN,warn)
%	FMIN	Minimize a function of one variable.
%	FMIN('F',x1,x2) attempts to return a value of x which is a local 
%	minimizer of F(x) in the interval x1 < x < x2.  'F' is a string 
%	containing the name of the objective function to be minimized.
%
%	FMIN('F',x1,x2,OPTIONS) uses a vector of control parameters.
%	If OPTIONS(1) is nonzero, intermediate steps in the solution are
%	displayed; the default is OPTIONS(1) = 0.  OPTIONS(2) is the termination
%	tolerance for x; the default is 1.e-4.  OPTIONS(14) is the maximum
%	number of steps; the default is OPTIONS(14) = 500.  The other components
%	of OPTIONS are not used as input control parameters by FMIN.  For more
%	information, see FOPTIONS.
%
%	FMIN('F',x1,x2,OPTIONS,P1,P2,...) provides for up to 10 additional
%	arguments which are passed to the objective function, F(X,P1,P2,...)
%
%	Example:
%	    fmin('cos',3,4) computes pi to a few decimal places.
%	    fmin('cos',3,4,[1,1.e-12]) displays the steps taken
%	    to compute pi to about 12 decimal places.  

%	Reference: "Computer Methods for Mathematical Computations",
%	Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.

%	Revised 10-14-88 JNL, 7-9-90 ACWG, 1-17-92, 6-15-92 CBM.
%	Original coding by Duane Hanselman, University of Maine.
% initialization

if nargin < 7
warn=0;
end

tol = options(2);
num = 1;
seps = sqrt(2.2204e-16);
c = 0.5*(3.0 - sqrt(5.0));
a = ax; b = bx;
v = a + c*(b-a);
w = v; xf = v;
d = 0.0; e = 0.0;
x= xf; fx = targetfn_iso(x,ratio,sigma,wN);  
fv = fx; fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;   tol2 = 2.0*tol1;

% Main loop

options(10) = 0; 
while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) )

    num = num+1;
    gs = 1;

    % Is a parabolic fit possible
    if abs(e) > tol1
        % Yes, so fit parabola
        gs = 0;
        r = (xf-w)*(fx-fv);
        q = (xf-v)*(fx-fw);
        p = (xf-v)*q-(xf-w)*r;
        q = 2.0*(q-r);
        if q > 0.0,  p = -p; end
        q = abs(q);
        r = e;  e = d;

        % Is the parabola acceptable
        if ( (abs(p)<abs(0.5*q*r)) & (p>q*(a-xf)) & (p<q*(b-xf)) )

            % Yes, parabolic interpolation step
            d = p/q;
            x = xf+d;
            step = '   num        xf        fx       parabolic';
     
            % f must not be evaluated too close to ax or bx
            if ((x-a) < tol2) | ((b-x) < tol2)
                si = sign(xm-xf) + ((xm-xf) == 0);
                d = tol1*si;
            end
        else
            % Not acceptable, must do a golden section step
            gs=1;
        end
    end
    if gs
        % A golden-section step is required
        if xf >= xm, e = a-xf;    else, e = b-xf;  end
            d = c*e;
            step = '   num        xf        fx       golden   ';
    end

    % The function must not be evaluated too close to xf
    si = sign(d) + (d == 0);
    x = xf + si * max( abs(d), tol1 );
    fu = targetfn_iso(x,ratio,sigma,wN);  

    % Update a, b, v, w, x, xm, tol1, tol2
    if fu <= fx
        if x >= xf, a = xf; else, b = xf; end
        v = w; fv = fw;
        w = xf; fw = fx;
        xf = x; fx = fu;
    else % fu > fx
        if x < xf, a = x; else,b = x; end
        if ( (fu <= fw) | (w == xf) )
            v = w; fv = fw;
            w = x; fw = fu;
        elseif ( (fu <= fv) | (v == xf) | (v == w) )
            v = x; fv = fu;
        end
    end
    xm = 0.5*(a+b);
    tol1 = seps*abs(xf) + tol/3.0; tol2 = 2.0*tol1;
    if num > options(14)
        if options(1)>=0
        disp('Warning: Maximum number of iterations has been exceeded');
        disp('       - increase options(14) for more iterations.')
        options(10)=num;
        return
        end
    end
end

fvm=fv;
options(8) = min(fv);
options(10)=num;
