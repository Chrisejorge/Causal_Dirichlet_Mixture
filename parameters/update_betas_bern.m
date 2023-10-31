function newBetas = update_betas_bern(Y,TX_matrix,z,k,beta_coeff, beta_coeff0,diagbetacov0)
% Update values of beta coefficient
% beta_coeffs must be n*d

% the number of y clusters
% k = length(unique(z(:,1)));

% Update values of beta_coeff for outcome regressions
newBetas = zeros(k,size(TX_matrix,2));
for jj = 1:k
    beta_coeff_current = beta_coeff(jj,:);
    tx_tmp = TX_matrix(z(:,1)==jj,:);
    y_tmp = Y(z(:,1)==jj);
    newBetas(jj,:) = update_beta_coeff_bern(beta_coeff0,diagbetacov0,beta_coeff_current,tx_tmp,y_tmp );
end %end for
% end
