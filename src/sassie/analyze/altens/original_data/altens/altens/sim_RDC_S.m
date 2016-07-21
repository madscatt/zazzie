function D=sim_RDC_S(vNH,S)
%----------------------------------------------------
%   df-nov-2011
%   simulate RDC for a given set of vectors and 
%   the alignment tensor S (not as diag and Euler separately)
%
%----------------------------------------------------
nres=size(vNH,1);
%--------- calculate matrix of direction cosines --------------
cosmatrix=vect2cosmat(vNH(:,2:4));

%--------- convert S into vector Q --------------------
Q=altens2vect(S);

%-------- generate RDCs --------------
D=cosmatrix*Q;
D=[vNH(:,1),D];

%----- input noise over Di ---------
%noise=ones(nres,1);
%noise=noise*0.01.*D(:,2)
%rand_noise=(rand(nres,1)*2)-1
%gen_noise=noise.*rand_noise
%D=[D,gen_noise]
%D=[D(:,1) D(:,2)+D(:,3)]
