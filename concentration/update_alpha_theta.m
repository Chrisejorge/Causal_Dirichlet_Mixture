function alpha_theta_new = update_alpha_theta(alpha_shape,alpha_rate,nobs,num_of_clusters, alpha_current)
% Define function to update alpha_theta prior
n =    nobs;
k = num_of_clusters;
eta = betarnd( alpha_current + 1, n, 1);

pieta = (alpha_shape + k - 1) / (alpha_shape + k - 1 + n * (1 - log(eta)) );

whichmix = binornd(1, pieta,1);
alpha_theta_new = whichmix * gamrnd( (alpha_shape + k), 1/(alpha_rate - log(eta)),1) ...
    + (1 - whichmix) * gamrnd((alpha_shape + k - 1), 1/(alpha_rate - log(eta)),1);
% end