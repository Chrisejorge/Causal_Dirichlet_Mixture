function newvar = update_covs_variance(nu0,tau0,c0,mu0,x_tmp)
% Define function to update variance parameter for covariate, N-Inv-chi^2 model
% x_tmp must be stored by column vectors

tmplen = size(x_tmp,1);
df_new = nu0 + tmplen;

if tmplen == 1
    varx = 0;
end

if tmplen > 1
    varx = var(x_tmp);
end

numer = nu0 * tau0 * tau0 + (tmplen - 1) .* varx + (c0 * tmplen/(c0 ...
    + tmplen)) * (mean(x_tmp) - mu0).^2;

newvar = sinvchi2rand(df_new, numer./df_new);
% end