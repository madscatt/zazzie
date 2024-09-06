function [vNH,R,IT]=lab2iner(filename,kheavy,kNH,kord,offs),
%--------------------------------------------------------------------
%	df-dec-97	df-aug-97	df-jul-96
%	get coordinates (and at.weight) of ALL or 
%	heavy atoms (N,C,O,S) from the pdb-data
%	calculate the inertia tensor and
%	define the inertia-tensor coordinate frame;
%	determine coordinates of the NH vectors in the 
%	inertia-tensor coordiante frame 
%
%	CALL: 	[vNH,R,IT]=lab2iner(filename,kheavy,kNH,kord,offs)
%	INPUT: 	filename (string)
%		kheavy (optional)=1 to get the heavy atoms only
%				 =0 to get all atoms   (default)
%		kNH (optional)=1 get & rot. NH vectors (default) 
%			      =0 skip the NH vect. rotation 
%			      =-1 skip the NH vect. at all
%		kord(optional)=1 IT in ascend.order
%			      =-1   in descending order
%			      =0    no addit.ordering (default)
%		offs	=0 if amide proton=HN (default)
%			=-1 if amide proton=H (e.g. insight, etc)
%	OUTPUT:	vNH -- [nres,vNHx,vNHy,vNHz] NH vect.
%			  coords in the inert-tens.frame
%		R  -- rotation matrix, lab2inert, (3x3)
%		IT -- inertia tensor  (3x3)
%
%	Note: use offs=-1 if reading a pdb-file with amide proton=H
%	uses:   atnamlst.m, readpdb.m (readasci.m),select.m,get_atwt.m 
%		inertens.m (row2mtrx.m,ordmtrx.m,sortmat.m),NHvect.m
%---------------------------------------------------------------------
format compact
%------- check input---------
if nargin<5, offs=0; 	end		%default: amide proton=HN
if nargin<4, kord=0; 	end		%default
if nargin<3, kNH=1; 	end		%default
if nargin<2, kheavy=0; 	end		%default
if kord~=1&kord~=0&kord~=-1, error('wrong kord!'); end
reslst=[];				%take all residues
%-------make atom-names list-------
if kheavy==1, atlst='hv'; else, atlst=[]; end 
[atlst,atdim]=atnamlst(atlst);
%-------read in the data set---------
[coor,atnam,at_res]=readpdb(filename,reslst,atlst,1);
%-------select atoms-----------------
if kheavy==1,
   sel=select(atnam,at_res,reslst,atlst,atdim);
else
   sel=[1:size(coor,1)];
end
nat=length(sel);
if nat==0, error('no atoms selected!!! wrong res/at list'); end
disp([num2str(nat),' atoms selected']);
%---------calc. atom Mw------------
atwt=get_atwt(atnam(sel,:));
disp('protein weight')
weight=sum(atwt)
%---------calc. Inertia Tensor-----
[R,IT,Rord]=inertens(coor(sel,:),atwt,-1);
ITdiag=diag(IT);
disp(['princ.values of inertia tensor: ',num2str(ITdiag(1)),', ',...
	num2str(ITdiag(2)),', ',num2str(ITdiag(3))]);
%------------get N,HN coord., calc.NH-vector---------------
if kNH~=-1,
  vNH=NHvect(coor,atnam,at_res,reslst,offs);
  %-----------NH-vect.coord. in the inert.-frame-----------
  if kNH==1,					%rotate NH
	vNH(:,2:4)=(Rord*R'*vNH(:,2:4)')';   %here I took into account that
	%the eig-procedure gives a rotation mtrx V comprised of 
	%eigen-vector columns instead or rows, i.e. V is supposed
	%to act on a row: X*V, ==> to rotate a column we need V',
	%hence R' instead of R
  else
     disp('the NH vector has NOT been rotated!!!');
  end
else
  vNH=[];
end
return
%========================================================
