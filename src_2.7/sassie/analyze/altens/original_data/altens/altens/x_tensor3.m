function [rot,S_val,S]=x_tensor3(a_inv,d_eff)

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

	v=real(v);
    
    d=diag(real(d));
  
    [dd i]=sort(abs(d));
   
    d=d(i);
  
 
    v=v(:,i);
  
    
    
    rot=v;
    
   
   	S_val=[d(1) d(2) d(3)];
	

return

%====================================================
