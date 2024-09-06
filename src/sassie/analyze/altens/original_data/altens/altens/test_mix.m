%test RDC mixing
S12=[9.3104 10.1584];
Euler=[0 90 0];
D1=sim_RDC(vNH,S12,Euler);
V=rotation_matrix(rot(1),rot(2),rot(3))';
D2=sim_RDC([vNH(:,1),vNH(:,2:4)*V],S12,Euler);
subplot(211)
bar([D1(:,1);D1(:,1)+50],[D1(:,2);D2(:,2)])
subplot(212)
bar([D1(:,1);D1(:,1)+50],[D1(:,2);D1(:,2)+D2(:,2)])