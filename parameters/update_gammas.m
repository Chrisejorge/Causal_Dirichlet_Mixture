function newGammas = update_gammas(T,X_matrix,z,k,kjs,gamma_coeff,gamma_coeff0,diaggammacov0)
% update values of gamma coefficient

% the number of y clusters
% k = length(unique(z(:,1)));

% kjs = zeros(k,1);
% for jj = 1:k
%     kjs(jj) = length(unique(z(z(:,1) == jj,2)));
% end

% numTclusters = size(unique(z(:,1:2),'rows'),1);
% newGammas = zeros(numTclusters,size(X_matrix,2));

newGammas = NaN(size(gamma_coeff));
count = 1;
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        gamma_coef_current = gamma_coeff(count,:);
        tmpidx = z(:,1) == jj & z(:,2) == l;
        x_tmp = X_matrix(tmpidx,:);
        t_tmp = T(tmpidx,:);
        newGammas(count,:) = update_gamma_coeff(gamma_coeff0,diaggammacov0,gamma_coef_current,x_tmp,t_tmp);
        count = count+1;
    end
end %end for

% end
