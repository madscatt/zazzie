%makedemo
freq=600.13;
TAUc=8.0;
Dz2Dx=1.5;
Dy2Dx=1.03;
AL=70;
BE=60;
GA=170;
th=rand(50,1)*180;
ph=rand(50,1)*360;
vNH=[];
vNH=[[1:50]',cos(ph*pi/180).*sin(th*pi/180)];
vNH=[vNH,sin(ph*pi/180).*sin(th*pi/180)];
vNH=[vNH,cos(th*pi/180)];
cd ../dynamics
R1=[];R2=[];R3=[];
for ii=1:50,
[r1,r2,r3]=reldata(freq,TAUc,1,0,0.02,-1,[Dz2Dx AL BE GA Dy2Dx],vNH(ii,2:4));
R1=[R1;[ii,r1(2:3)]];
R2=[R2;[ii,r2(2:3)]];
R3=[R3;[ii,r3(2:3)]];
end
%--- add noise
R1(:,2)=R1(:,2).*(1+randn(50,1)*0.02);
R2(:,2)=R2(:,2).*(1+randn(50,1)*0.02);
R3(:,2)=R3(:,2).*(1+randn(50,1)*0.02);
cd ../RotDif;
rotdif(freq,R1,R2,R3,vNH,[],'demo');