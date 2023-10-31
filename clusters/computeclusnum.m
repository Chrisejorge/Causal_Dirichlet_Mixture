function [k,kjs,kljs] = computeclusnum(z)
% Calculation of the number of clusters

% Calculate number of y-clusters
k = length(unique(z(:,1)));
% Calculate number of t-cluster within each y-cluster
kjs = zeros(k,1);
for jj = 1:k
    kjs(jj) = length(unique(z(z(:,1) == jj,2)));
end

% Calculate number of x-cluster within each y-cluster, t-cluster
maxNumTclusters = max(z(:,2));
kljs = zeros(k,maxNumTclusters);
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        tmpidx = z(:,1) == jj & z(:,2) == l;
        klj = length(unique(z(tmpidx,3)));
        kljs(jj,l) = klj;
    end
end %end for

% end