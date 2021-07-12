function vNH=getNHvect(coor,atnam,at_res,reslst,offs)
%-------------------------------------------------
%  df-dec-97
%       retrieve the NH vectors from coor data set
%       offs=-1 if amide hydr.=H in the pdb-data set
%       uses: select.m
%-------------------------------------------------
  if nargin<5, offs=0;  end             %default
  %----------select the amide N's-------
  selN=select(atnam,at_res,reslst,'N ',2);
  iN=length(selN);
  %----------select the amide H's-------
  if offs==-1,                          %insight-type naming
    selH=select(atnam,at_res,reslst,' H ',3,-1);
  else
    selH=select(atnam,at_res,reslst,'HN',2,offs);
  end
  iH=length(selH);
  if iH==0, 
     disp(['offs=',num2str(offs),': no amide protons found!']);
     vNH=[];
     return   
  end
  disp([num2str(iH),' amide protons found']);
  vNH=NaN*ones(iH,4);
  vNH(:,1)=at_res(selH,2);
  %-----filter the PROlines out ------
  for ii=1:iH,
    for jj=1:iN,
       if at_res(selN(jj),2)==at_res(selH(ii),2),
          vNH(ii,2:4)=coor(selH(ii),:)-coor(selN(jj),:);
          break
       end
    end
  end
  %------normalize the NH vectors------
  for ii=1:iH,
    vNH(ii,2:4)=vNH(ii,2:4)/(sqrt(vNH(ii,2:4)*vNH(ii,2:4)'));
  end
return
%=====================================================

