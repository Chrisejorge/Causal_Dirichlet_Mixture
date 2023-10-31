function [njs,nljs,nhljs] = computenuminclus(z,k,kjs,kljs)
% Calculation of numbers of data in each cluster

maxNumTclusters = max(z(:,2));
njs = zeros(k,1);
for jj = 1:k
    nj = size(z(z(:,1) == jj),1);
    njs(jj) = nj;
end

maxNumXclusters = max(z(:, 3));
nljs =  zeros(k,maxNumTclusters);
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        tmpidx = z(:,1) == jj & z(:,2) == l;
        nljs(jj,l) = size(z(tmpidx),1);
    end
end %end for jj

nhljs = zeros(k,maxNumTclusters,maxNumXclusters);
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        klj = kljs(jj,l);
        for h = 1:klj
            tmpidx = z(:,1) == jj & z(:,2) == l & z(:,3) == h;
            nhljs(jj,l,h) = size(z(tmpidx),1);
        end
    end
end %end for jj
% end