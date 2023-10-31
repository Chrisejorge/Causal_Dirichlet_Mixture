function clusValue = init_cluster_hclus(Y,TX,X,ky,kt,kx)
% Initialize cluster memberships by Hierarchical Clustering
% ky y-clusters
% kt t-clusters within each y-cluster
% kx x-clusters within each t-cluster
% Set initial values for cluster membership,y,t,x

n = size(Y,1);
% Use k-means to find 3 whole data clusters
clusValue = ones(n,3);

% Init  y-cluster
clusValue(:,1) = clusterdata([TX,Y],'Linkage','ward','Maxclust',ky);

% Make kt t-clusters within each y-cluster
for jj = 1:ky
    if sum(clusValue(:,1) == jj) > kt 
           clusValue(clusValue(:,1) == jj, 2) = clusterdata(TX(clusValue(:,1) == jj,:),'Linkage','ward','Maxclust',kt);
    else
            clusValue(clusValue(:,1) == jj, 2) = clusterdata(TX(clusValue(:,1) == jj,:),'Linkage','ward','Maxclust',sum(clusValue(:,1) == jj));
    end
end % end for

% Make kx x-clusters within each t-cluster
for jj = 1:ky
    for l = 1:kt
        tmpidx = clusValue(:,1) == jj & clusValue(:,2) == l;
        if  sum(tmpidx) > kx
            clusValue(tmpidx, 3) = clusterdata(X(tmpidx,:),'Linkage','ward','Maxclust',kx);
        else
             clusValue(tmpidx, 3) = clusterdata(X(tmpidx,:),'Linkage','ward','Maxclust',sum(tmpidx));
        end
    end % end for
end
