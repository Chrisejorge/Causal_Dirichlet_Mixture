function newmean = update_covs_mean(c0,mu0,x_tmp, tau2_current)
% Define function to update mean parameter for covariate
% x_tmp must be stored by column vectors
% tau2_current must be stored by row vector

tmplen = size(x_tmp,1);
newvar = 1./(c0./tau2_current + tmplen./tau2_current);
%newmean = (mu0 * c0/tau2_current + mean(x_tmp) * length(x_tmp)/tau2_current) * newvar;
newmean = (mu0 * c0 + mean(x_tmp) .* tmplen) / (c0 + tmplen);
newmean = normrnd(newmean, sqrt(newvar));
% end