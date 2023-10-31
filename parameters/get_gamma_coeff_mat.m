function gamma_coeff_mat = get_gamma_coeff_mat(z,k,kjs,gamma_coeff)
% Get matrix for gamma

maxNumTclusters = max(z(:,2));
gamma_coeff_mat = zeros(k+1,maxNumTclusters+1,size(gamma_coeff,2));
count = 1;
for jj = 1:k
    kj = kjs(jj);
    for l = 1:kj
        gamma_coeff_mat(jj,l,:) = gamma_coeff(count,:);
        count = count+1;
    end
end % jj

% end