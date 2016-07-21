altens([2:7,12:17,22:32,41:45,48,49,58:59,63:70],vNH_opD_d,RDC(1:76,:));

SS=[8.9875   18.3275];

RDC_opD_d=sim_RDC(vNH_opD_d,SS,[287.849 89.917 64.982]);
RDC_opD_p=sim_RDC(vNH_opD_p,SS,[287.849 89.917 64.982]);
bar([RDC_opD_d(:,1);RDC_opD_p(:,1)+100],[RDC_opD_d(:,2);RDC_opD_p(:,2)])
RDC_1aar_d=sim_RDC(vNH_1aar_d,SS,[287.849 89.917 64.982]);
RDC_1aar_p=sim_RDC(vNH_1aar_p,SS,[287.849 89.917 64.982]);

bar([RDC_1aar_d(:,1);RDC_1aar_p(:,1)+100],[RDC_1aar_d(:,2);RDC_1aar_p(:,2)])
cd ../df_helpers/
RDC_d=combine2(RDC_1aar_d,RDC_opD_d,2);
edit combine2
RDC_d=combine2(RDC_1aar_d,RDC_opD_d,RDC_opD_d(:,1),2);
RDC_p=combine2(RDC_1aar_p,RDC_opD_p,RDC_opD_p(:,1),2);
figure
bar([RDC_d(:,1);RDC_p(:,1)+100],[RDC_d(:,2)+RDC_d(:,3);RDC_p(:,2)]+RDC_p(:,3)])
bar([RDC_d(:,1);RDC_p(:,1)+100],[RDC_d(:,2)+RDC_d(:,3);RDC_p(:,2)+RDC_p(:,3)])
bar([RDC_d(:,1);RDC_p(:,1)+100],[RDC_d(:,2)+RDC_d(:,3);RDC_p(:,2)+RDC_p(:,3)]/2)
p=0.5;bar([RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)]/2)
p=0.5;bar([RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)])
p=0.7;bar([RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)])

p=0.6;bar([RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)])

p=0.6;RDC_pp=[[RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)]];
altens([2:7,12:17,22:32,41:45,48,49,58:59,63:70],vNH_opD_d,RDC_pp(1:61,:));
altens([2:7,12:17,22:32,41:45,48,49,58:59,63:70]+100,[vNH_opD_p(:,1)+100,vNH_opD_p(:,2:4)],RDC_pp(62:end,:));
altens([[2:7,12:17,22:32,41:45,48,49,58:59,63:70],[2:7,12:17,22:32,41:45,48,49,58:59,63:70]+100],[vNH_opD_d;[vNH_opD_p(:,1)+100,vNH_opD_p(:,2:4)]],RDC_pp);

p=0;RDC_pp=[[RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)]];


p=0.3;RDC_pp=[[RDC_d(:,1);RDC_p(:,1)+100],[p*RDC_d(:,2)+(1-p)*RDC_d(:,3);p*RDC_p(:,2)+(1-p)*RDC_p(:,3)]];
>> RRR=combine2(RDC_pp,RDC,RDC_pp(:,1),2);
>> plot(RRR(:,2),RRR(:,3),'o')