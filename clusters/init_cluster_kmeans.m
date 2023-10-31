function clusValue = init_cluster_kmeans(Y,TX,X,ky,kt,kx)
% Initialize cluster memberships by kmeans
% ky y-clusters
% kt t-clusters within each y-cluster
% kx x-clusters within each t-cluster
% Set initial values for cluster membership,y,t,x

n = size(Y,1);
% Use k-means to find 3 whole data clusters
clusValue = ones(n,3);

% Init  y-cluster
clusValue(:,1) = kmeans([TX,Y],ky);

% Make kt t-clusters within each y-cluster
for jj = 1:ky
    if sum(clusValue(:,1) == jj) > kt 
        clusValue(clusValue(:,1) == jj, 2) = kmeans(TX(clusValue(:,1) == jj,:),kt);
    else
        clusValue(clusValue(:,1) == jj, 2) = kmeans(TX(clusValue(:,1) == jj,:),sum(clusValue(:,1) == jj));
    end
end % end for

% Make kx x-clusters within each t-cluster
for jj = 1:ky
    for l = 1:kt
        tmpidx = clusValue(:,1) == jj & clusValue(:,2) == l;
        if  sum(tmpidx) > kx
            clusValue(tmpidx, 3) = kmeans(X(tmpidx,:),kx);
        else
            clusValue(tmpidx, 3) = kmeans(X(tmpidx,:),sum(tmpidx));
        end
    end % end for
end

% end

