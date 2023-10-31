function [gamma_coeff0,diaggammacov0,newGammas] = init_gamma_params(T,X,X_matrix,z,k,kjs)
% initialize parameters of gamma coefficients

% logit regression,T~X
gamma_coeff0 = (glmfit(X,T,'binomial'))';

diaggammacov0 = diag(4*ones(size(X_matrix,2),1));

numTclusters = size(unique(z(:,1:2),'rows'),1);
gamma_coeffs_current = NaN(numTclusters,size(X_matrix,2));

count = 0;
for jj = 1:k
    tmp = mvnrnd(gamma_coeff0,diaggammacov0,1);
    kj = kjs(jj);
    gamma_coeffs_current(count+1:count+kj,:) = repmat(tmp,kj,1);
    count = count + kj;
end %end for

newGammas = update_gammas(T,X_matrix,z,k,kjs,gamma_coeffs_current,gamma_coeff0,diaggammacov0);

% end