%this script takes RDC data for both distal and proximal domains, 
%rotates the proximal domain on a grid in Euler angles {alpha,beta,gamma},
%and for each point runt ALTENS and records the chi^2 value
%vNH=coordinates for the distal domain; 
%vNHp=starting coordinates for the proximal domain (make res# + 100 to distinguish from the distal!)
%don't forget to define core residues!

al=[0:1:180];
be=[0:1:180];  %be=[0:5:180];
ga=[0:1:180];  %ga=[0:5:180];


%al=[0:1:359];
%be=[15:1:25];
%ga=[70:1:90];

%al=[170:.1:190];

%al=[0:1:10];
%be=[176:1:180];
%ga=[120:1:150];

%core=[2:9,11:71,102:171];

%optimized core: D: [2:7,9,11:30,32:34,36:71]; P: [2:6,8:11,13:22,24:49,51:71]
%core=[coreD_rdc,coreP_rdc];
%core=[coreD,coreP];
%RDC=[rdcD;rdcPs];
%RDC=[rdcD;rdcP];
RR=[];
RR=zeros(length(ga)*length(be)*length(al),4);
ind=0;
indt=0;
for ii=1:length(ga), 
    %temp=[];
    temp=zeros(length(be)*length(al),1);
    for jj=1:length(be), 
        for kk=1:length(al),
           [a,b,R]=altens_fly(core,[vNH;[vNHp(:,1),rotate_vectors(vNHp(:,2:4),[al(kk),be(jj),ga(ii)],[0 0 0],0)]],RDC);
           ind=ind+1; 
           indt=indt+1; 
           %RR=[RR;[al(kk),be(jj),al(ii),R]];
           RR(ind,:)=[al(kk),be(jj),ga(ii),R];
           %temp=[temp,R]; 
	     temp(indt)=R; 
        end
    end
    disp(min(temp))
    indt=0;
end
