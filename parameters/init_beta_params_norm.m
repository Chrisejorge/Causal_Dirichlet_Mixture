function [beta_coeff0,diagbetacov0,precbetacov0,sigma2_coeff,beta_coeff] = init_beta_params_norm(Y,TX,TX_matrix,z,k,sigma2_shape,sigma2_scale)
% Update values of parameters
% for outcome regression coefficients
% Parameters in outcome regression
% linear regression,Y~T,X
beta_coeff0 = (glmfit(TX,Y,'normal'))';

diagbetacov0 = diag(4*ones(size(TX_matrix,2),1));

precbetacov0 = inv(diagbetacov0);

[sigma2_coeff, beta_coeff ] = update_sigma2s_betas(Y,TX_matrix,z,k,beta_coeff0,precbetacov0,sigma2_shape,sigma2_scale);

