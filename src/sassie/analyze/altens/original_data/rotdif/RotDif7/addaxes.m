function [coor,atnam,at_res]=addaxes(fname,fname_output,rot_angle,rad,reslst,atlst,model)
%-----------------------------------------------------------------
%   df-apr-13 if rot_angle has only 2 entries, treat as axial model
%   df-march-2002
%   modified by ow based on the program readpdb of david fushman
%	read and transform pdb datra set according to a given rotation matrix
%   uses: readasci.m
%  fname : name of the pdb input file
%  fname_output : name of the pdb output file
%  rot_angle : vector containing the three euler angles
%
%Copyright (c) 2004 by David Fushman, University of Maryland
%-----------------------------------------------------------------
if nargin<7, model=1; 	end		%default: first model
if nargin<6, atlst=[];  end     %default: all atoms
if nargin<5, reslst=[]; end		%default: all resid.
kax=0;
if length(rot_angle)==2, kax=1; rot_angle=[rot_angle,0]; end  %draw only one axis
pdb=readasci(fname);
disp(['got the ',fname,' data set, analyzing...']);
nlin=size(pdb,1);
select_at=zeros(nlin,1);
header = [];
isel=1;
nmod=1;
termflag=0;				%term flag off
for ii=1:nlin,
  if strcmp('TER',pdb(ii,1:3))|strcmp('END',pdb(ii,1:3))|strcmp('.',pdb(ii,1)), 
     if termflag==0,
        nmod=nmod+1;
        termflag=1;			%term flag on
        if nmod>model, break; end
     end
  elseif strcmp('HEADER',pdb(ii,1:6)), 
      header = [header; pdb(ii,:)]; 
  else
     termflag=0;			%term flag off
     if nmod==model			%read in the structure
       if strcmp('ATOM',pdb(ii,1:4)), 
         if ~isempty(reslst),		%select res#
            if ~isempty(find(reslst==str2num(pdb(ii,23:26)))),
               select_at(isel)=ii; isel=isel+1;
            end
         else  select_at(isel)=ii; isel=isel+1;
         end
       end
     end
  end
end
sel=select_at(1:isel-1);
nat=length(sel);
if nat==0, error('no atoms found!!! wrong filename or model'); end  
disp([num2str(nat),' atoms read in']);
atnam=pdb(sel,1:26);
coor=zeros(nat,3);
at_res=zeros(nat,2);
at_res(:,1)=str2num(pdb(sel,7:11));
at_res(:,2)=str2num(pdb(sel,23:26));

fprintf('rotating axis with alpha = %6.3f ; beta = %6.3f ; gamma = %6.3f \n',rot_angle(1),rot_angle(2),rot_angle(3));

%----- extract coordinates of pdb file -----------------
coor(:,1:3)=str2num(pdb(sel,32:54));


%----- center of mass -----------
gc=mean(coor);

%----- extract the three euler angles ------------------
alpha=rot_angle(1)*pi/180.0;
beta=rot_angle(2)*pi/180.0;
gamma=rot_angle(3)*pi/180.0;




%----- create rotation matrix and rotate coordinate-----
rot=rotation_matrix(alpha,beta,gamma);

vect1=rot(1,:).*rad;
vect2=rot(2,:).*rad;
vect3=rot(3,:).*rad;

ca=gc;				%CA-- in the geometrical center

vx=(gc+vect1);
vy=(gc+vect2);
vz=(gc+vect3);
vmz=(gc-vect3);


%----- extract end of pdb file -------------------------
fre=zeros(nat,2);
fre(:,1:2)=str2num(pdb(sel,55:66));
ter=pdb(sel,67:end);
terwidth = size(ter,2);
if terwidth >= 14,
    ter = ter(:,1:14);
end


%----- print new pdb file in output file ---------------
fid=fopen(fname_output,'w');


%----- header of pdb file (may be modified) ------------
addheader = ['HEADER    MODIFIED: DIFFUSION AXES ADDED          ',upper(date),'            ADDAXES\n'];
fprintf(fid,addheader);
if ~isempty(header),
    for ii=1:size(header,1),
      headertmp = deblank(header(ii,:));
      headertmp = [headertmp,blanks(14-length(headertmp)),'\n'];
      fprintf(fid,headertmp);
    end
end    

%----- new coordinates of atoms ------------------------
for jj=1:nat
    %padding with blanks
    tertmp = deblank(ter(jj,:));
    tertmp = [tertmp,blanks(14-length(tertmp))];
    fprintf(fid,'%s%12.3f%8.3f%8.3f%6.2f%6.2f%s\n',atnam(jj,1:26),coor(jj,1),coor(jj,2),coor(jj,3),fre(jj,1),fre(jj,2),tertmp);
end


printline='  GLY   999    %8.3f%8.3f%8.3f  1.00 32.83      DIFFAXIS\n';
  fprintf(fid,['ATOM   9995  N ',printline],vz(1),vz(2),vz(3));
  fprintf(fid,['ATOM   9996  CA',printline],ca(1),ca(2),ca(3));
  fprintf(fid,['ATOM   9997  C ',printline],vmz(1),vmz(2),vmz(3));
  if kax==0,       %do not draw X,Y axes
    fprintf(fid,['ATOM   9998 1HA',printline],vx(1),vx(2),vx(3));
    fprintf(fid,['ATOM   9999 2HA',printline],vy(1),vy(2),vy(3));
  end
%----- end of pdb file ---------------------------------

fprintf(fid,'END                                                                      ADDAXES\n');

fclose(fid);






return
%=================================================================
