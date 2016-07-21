function D=sim_align_tens(NH,S)


% simulate RDC with a given set of vectors and S parameters



nres=size(NH,1);

%--------- calculate direction cosine --------------

z=dir_cos(NH(:,2:4));


%-------- generate matrix Di --------------

D=gen_Di(z,S);
D=[NH(:,1),D]


%----- input noise over Di ---------
noise=ones(nres,1);
noise=noise*0.01.*D(:,2)


rand_noise=(rand(nres,1)*2)-1

gen_noise=noise.*rand_noise

D=[D,gen_noise]

D=[D(:,1) D(:,2)+D(:,3)]



