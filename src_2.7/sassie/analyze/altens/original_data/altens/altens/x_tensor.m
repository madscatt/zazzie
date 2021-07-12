function [rot,S_val,S]=x_tensor(a_inv,d_eff)

q_v=a_inv*d_eff;

%--------Reconstruct the Q-matrix---------------
   
   S(1,1)=-q_v(1)-q_v(2);
   S(2,2)=q_v(1);
   S(3,3)=q_v(2);
   S(1,2)=q_v(3);
   S(1,3)=q_v(4);
   S(2,3)=q_v(5);
   S(3,1)=S(1,3);
   S(3,2)=S(2,3);
   S(2,1)=S(1,2);
   
%---------Diagonalize the Q matrix--------------

[v,d]=eig(S);


%---------Sort the eigensystem------------------

	Szz=+inf;
	Syy=0.0;
	Sxx=-inf;
	for i=1:3
   	if d(i,i) < Szz
      	Szz=d(i,i);
      	rot(:,3)=v(:,i);
   	end
   	if d(i,i) > Sxx
      	Sxx=d(i,i);
      	rot(:,1)=v(:,i);
   	end
	end
	for i=1:3
   	if d(i,i) ~= Sxx & d(i,i) ~= Szz
      	Syy=d(i,i);
      	rot(:,2)=v(:,i);
   	end
end
   S_val=[Sxx Syy Szz];
	

return

%====================================================
