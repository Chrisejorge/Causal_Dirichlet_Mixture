function gamma_coeff_draw = update_gamma_coeff(gamma_coeff0,diaggammacov0, gamma_coef_current,x_tmp,t_tmp)
% Define function to update gamma coefficients in regression.

tmpcov0 = diag(0.01*ones(size(gamma_coeff0,2),1));

proposed = mvnrnd(gamma_coef_current,tmpcov0,1);

loglike_prop = sum( log(binopdf(t_tmp,1,expit(x_tmp*proposed')) + eps) )+...
    logmvnpdf(proposed,gamma_coeff0,diaggammacov0) + ...
    logmvnpdf(gamma_coef_current,proposed,tmpcov0);

loglike_old = sum( log(binopdf(t_tmp,1,expit(x_tmp*gamma_coef_current')) + eps) ) ...
    +logmvnpdf(gamma_coef_current,gamma_coeff0,diaggammacov0) ...
    + logmvnpdf(proposed,gamma_coef_current,tmpcov0);


tmp = min(exp(loglike_prop-loglike_old),1);

% fprintf(['exp ' num2str(exp(loglike_prop-loglike_old))   '\n']);

if rand < tmp
    gamma_coeff_draw = proposed;
    %          fprintf(['gamma proposed ' num2str(loglike_prop)   '\n']);
    %fprintf(['gamma exp ' num2str(exp(loglike_prop-loglike_old))   '\n']);
else
    %          fprintf([' gamma gamma_coef_current ' num2str(gamma_coef_current)   '\n']);
    gamma_coeff_draw = gamma_coef_current;
end
% end
