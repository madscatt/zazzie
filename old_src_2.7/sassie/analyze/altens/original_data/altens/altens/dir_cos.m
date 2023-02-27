function z=dir_cos(coordinates)


ncoor=size(coordinates,1);

coordinates=coordinates';
z=[];
%------ define X, Y and Z axis ------------

X=[1 0 0];
Y=[0 1 0];
Z=[0 0 1];

X=X';Y=Y';Z=Z';


%------ calculate angles with respect to each axis (in radians) -----

%for ii=1:ncoor
 %   z(ii,1)=subspace(coordinates(:,ii),X);
  %  z(ii,2)=subspace(coordinates(:,ii),Y);
   % z(ii,3)=subspace(coordinates(:,ii),Z);
    
    
   %end

for ii=1:ncoor
    z(ii,1)=acos(coordinates(1,ii));
    z(ii,2)=acos(coordinates(2,ii));
    z(ii,3)=acos(coordinates(3,ii));
end

%for ii=1:ncoor
 %   z(ii,1)=atan(norm(cross(coordinates(:,ii),X))/dot(coordinates(:,ii),X));
  %  z(ii,2)=atan(norm(cross(coordinates(:,ii),Y))/dot(coordinates(:,ii),Y));
   % z(ii,3)=atan(norm(cross(coordinates(:,ii),Z))/dot(coordinates(:,ii),Z));
    
  %  z=[z;z(ii,:)];
  %end


%z


return