function [newSigma2s, newBetas ] = update_sigma2s_betas(Y,TX_matrix,z,k,beta_coeff0,precbetacov0,sigma2_shape,sigma2_scale)
% update values of beta coefficient and sigma2 coefficient for outcome regressions

newSigma2s = zeros(k,1);
newBetas = zeros(k,size(TX_matrix,2));
for jj = 1:k
    tmpidx = jj;
    tx_tmp = TX_matrix(z(:,1)==tmpidx,:);
    y_tmp = Y(z(:,1)==tmpidx);
    precbetacov = precbetacov0' + tx_tmp' * tx_tmp; % inv(A)' = (inv(A))'
    betaLambda = inv(precbetacov);
    betaU = precbetacov \ (precbetacov0'*beta_coeff0' + tx_tmp' * y_tmp);
    newBetas(jj,:) = mvnrnd(betaU',betaLambda',1);
    an = sigma2_shape + length(y_tmp)/2;
    bn = sigma2_scale + 0.5 * (y_tmp'*y_tmp + beta_coeff0 * precbetacov0' * beta_coeff0'-betaU'*precbetacov*betaU);
    newSigma2s(jj,:) = 1.0/gamrnd(an,1/bn);
end
