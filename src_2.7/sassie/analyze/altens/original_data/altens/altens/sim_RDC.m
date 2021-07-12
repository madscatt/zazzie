function D=sim_RDC(vNH,Sxxyy,Euler)
%----------------------------------------------------
%   df-nov-2011
%   simulate RDC for a given set of vectors and 
%   alignment tensor: Sxxyy = [Sxx,Syy] diag elements
%   Euler = {alpha, beta, gamma} in degrees
%
%----------------------------------------------------
nres=size(vNH,1);
%--------- calculate matrix of direction cosines --------------
cosmatrix=vect2cosmat(vNH(:,2:4));

%---------- reconstruct alignment tensor S ----------------
Sdiag=diag([Sxxyy,-sum(Sxxyy)]);
alpha = Euler(1)*pi/180;
beta = Euler(2)*pi/180;
gamma = Euler(3)*pi/180;

V = rotation_matrix(alpha,beta,gamma)';
S = V*(Sdiag*transpose(V))


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
