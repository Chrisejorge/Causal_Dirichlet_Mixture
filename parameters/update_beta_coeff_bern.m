function beta_coeff_draw = update_beta_coeff_bern(beta_coeff0,diagbetacov0,beta_coeff_current,tx_tmp,y_tmp)
% Define function to update beta coefficients in regression. Input current value
% of beta coefficient, a vector of predictors(including intercept) and outcome.
% Currently uses Normal(mean=0,sd=2) prior for beta coefficients. Change this
% within the function if needed.
% metropolis hastings acceptence
% p(beta|y,x) \propto p(y|beta,x) p(beta)
% p(beta|y,x) p(beta^{*}|beta)
% p(y|beta,x) = \prod p(y_{i}|x_{i},beta)
tmpcov0 = diag(0.01*ones(size(beta_coeff0,2),1));

proposed = mvnrnd(beta_coeff_current,tmpcov0,1);

loglike_prop = sum( log(binopdf(y_tmp,1,expit(tx_tmp*proposed')) + eps ) )+...
    logmvnpdf(proposed,beta_coeff0,diagbetacov0) + ...
    logmvnpdf(beta_coeff_current,proposed,tmpcov0);

loglike_old = sum( log(binopdf(y_tmp,1,expit(tx_tmp*beta_coeff_current')) + eps) ) ...
    +logmvnpdf(beta_coeff_current,beta_coeff0,diagbetacov0) ...
    + logmvnpdf(proposed,beta_coeff_current,tmpcov0);

tmp = min(exp(loglike_prop-loglike_old),1);

if rand(1) < tmp
    beta_coeff_draw = proposed;
    %         fprintf(['beta proposed ' num2str(loglike_prop)   '\n']);
else
    %          fprintf(['beta beta_coeff_current ' num2str(beta_coeff_current)   '\n']);
    beta_coeff_draw = beta_coeff_current;
end
% end