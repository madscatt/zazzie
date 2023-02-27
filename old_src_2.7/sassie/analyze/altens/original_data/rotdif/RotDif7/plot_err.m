function plot_err(X,Y,matrix,horiz,vert)

% Plot the surface representing errors
% X: absissa, Y: ordinate, matrix: errors defining surface
% horiz: a vector containing the values of sigma to be projected
% vert: a vector containing the values of Dyy/Dxx to be projected



%---- plot surface -----

figure
surf(Y,X,matrix)
shading interp
alpha 0.7


ind_horiz=[];
ind_vert=[];

%--- find index ----
for ii=1:length(horiz),
ind_hor=find(X(:,1)==horiz(ii));
ind_horiz=[ind_horiz,ind_hor];
end

for ii=1:length(vert),
ind_ver=find(Y(:,1)==vert(ii));
ind_vert=[ind_vert,ind_ver];
end


ind_vert
ind_horiz

l_horiz=length(ind_horiz);
l_vert=length(ind_vert);
if ~isempty(horiz),
for ii=1:l_horiz
    hold on 
    plot3(Y,0.1*ones(10,1),matrix(ind_horiz(ii),:),'k','LineWidth',2);
    hold on
    plot3(Y,ones(10,1)*X(ind_horiz(ii),1),matrix(ind_horiz(ii),:),'.k','MarkerSize',12);
    hold on
    for jj=1:10
        hold on
        plot3([Y(jj) Y(jj)],[X(ind_horiz(ii)) 0.1],[matrix(ind_horiz(ii),jj) matrix(ind_horiz(ii),jj)],'--k');
    end
end
end
    
if ~isempty(vert), 
for ii=1:l_vert
    hold on 
    plot3(1.2*ones(10,1),X,matrix(:,ind_vert(ii)),'k','LineWidth',2);
    hold on
    plot3(Y(ind_vert(ii),1)*ones(10,1),X,matrix(:,ind_vert(ii)),'.k','MarkerSize',12);
    hold on
    for jj=1:10
        hold on
        plot3([Y(ind_vert(ii)) 1.2],[X(jj) X(jj)],[matrix(jj,ind_vert(ii)) matrix(jj,ind_vert(ii))],'--k');
    end
end
end


