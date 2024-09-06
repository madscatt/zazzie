function [z, zsel] = euler2all(alpha,beta,gamma)
%-----------------------------------------------------
%   df-apr-04
%   for a given set of euler angles (in degrees), 
%   produce all other possible angles
%   based on the symmetry properties of
%   the diffusion tensor
%   {a,b,g}={a,b,g+180o}={a+180o,180o-b,180o-g}={a+180o,180o-b,360o-g}
%   and select the set that has all three angles in the [0 180]-cube
%------------------------------------------------------    
if nargin == 1, 
    if length(alpha) == 3,  %input Euler angles as a vector 
        beta = alpha(2); gamma = alpha(3); alpha = alpha(1);   
    else
        error('wrong input!');
    end   %default no selection
end
z = NaN*ones(4,3); zsel = NaN;
z(1,:) = [alpha,beta,gamma];
z(2,:) = [alpha,beta,gamma+180];
z(3,:) = [alpha+180,180-beta,180-gamma];
z(4,:) = [alpha+180,180-beta,360-gamma];

%convert angles to "normal" range
ind = find(z >= 270);
if ~isempty(ind), z(ind) = z(ind) - 360; end
ind = find(z <= -270);
if ~isempty(ind), z(ind) = z(ind) + 360; end

%select
indsel = find(z(:,1)>=0 & z(:,2)>=0 & z(:,3)>=0 & z(:,1)<=180 & z(:,2)<=180 & z(:,3)<=180);
if ~isempty(indsel), zsel = z(indsel,:); end

%======================================================
