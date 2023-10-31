function [beta_coeff0,diagbetacov0,newBetas] = init_beta_params_bern(Y,TX,TX_matrix,z,k)
% Initialize parameters of beta coefficients

% for outcome regression coefficients
% Parameters in outcome regression
% logit regression,Y~T,X
beta_coeff0 = (glmfit(TX,Y,'binomial'))';

diagbetacov0 = diag(4*ones(size(TX_matrix,2),1));

% Update values of beta_coeffs for outcome regressions
beta_coeffs_current = mvnrnd(beta_coeff0,diagbetacov0,k);

newBetas = update_betas_bern(Y,TX_matrix,z,k, beta_coeffs_current, beta_coeff0,diagbetacov0);

% end
