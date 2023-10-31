function [median_psi_rr,sd_psi_rr,cov_psi_rr,width_psi_rr] = compute_psi_rr_stats(psi_rr_true,psi_rr)
% Calculate parameters of psi_rr

% Calculate median of posterior for average relative risk
median_psi_rr = median(psi_rr(~isnan(psi_rr)));
sd_psi_rr = std(psi_rr(~isnan(psi_rr)));


% Determine if true relative risk is in credible interval
cov_psi_rr = (min(quantile(psi_rr(~isnan(psi_rr)),[0.025,0.975])) < psi_rr_true) & (psi_rr_true < max(quantile(psi_rr(~isnan(psi_rr)),[0.025,0.975])));

% Calculate width of credible interval
width_psi_rr = abs(max(quantile(psi_rr(~isnan(psi_rr)),[0.025,0.975]))-min(quantile(psi_rr(~isnan(psi_rr)),[0.025,0.975])) );
% end