function vectors1=rotate_vectors(vectors,rot_angles,frame_angles,kprint)
%-----------------------------------------------------------------
%   df-nov-2011
%   Rotate vectors by a given set of Euler angles (in degrees)
%   with respect to a frame rotated by by frame_angles
%   vectors = [x y z]
%-----------------------------------------------------------------
if nargin < 4, kprint = 1; end              %default: verbal
if nargin < 3, frame_angles=[0 0 0]; end    %default: rotation in the oritinal frame 

nat=size(vectors,1);
if kprint==1,
  fprintf('rotating vectors by alpha = %6.3f ; beta = %6.3f ; gamma = %6.3f \n',...
    rot_angles(1),rot_angles(2),rot_angles(3));
  fprintf('with respect to rotation frame defined by the Euler angles {%6.3f, %6.3f, %6.3f} \n',...
    frame_angles(1),frame_angles(2),frame_angles(3));
end
%----- prepare the euler angles ------------------
%---- for the rotation-----
alpha=rot_angles(1)*pi/180;
beta=rot_angles(2)*pi/180;
gamma=rot_angles(3)*pi/180;

%---- for the rotation frame-----
alphaFR=frame_angles(1)*pi/180.0;
betaFR=frame_angles(2)*pi/180.0;
gammaFR=frame_angles(3)*pi/180.0;

%----- create rotation matrices ------
Rot=rotation_matrix(alpha,beta,gamma);
RotFR=rotation_matrix(alphaFR,betaFR,gammaFR);

%express coordinates in the rotation frame, rotate, and then 
%express coordinates back in the x,y,z frame-----

%----- put domain along this orientation ------
vectors1=(RotFR'*(Rot)'*RotFR*(vectors'))'; 

return
%=================================================================
