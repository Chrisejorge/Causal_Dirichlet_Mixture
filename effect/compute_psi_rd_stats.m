function [median_psi_rd,sd_psi_rd,cov_psi_rd,width_psi_rd] = compute_psi_rd_stats(psi_rd_true,psi_rd)
% Calculate parameters of psi_rd

% Calculate median of posterior for average risk difference
median_psi_rd = median(psi_rd(~isnan(psi_rd)));
sd_psi_rd = std(psi_rd(~isnan(psi_rd)));

% Determine if true risk difference is in credible interval
cov_psi_rd = (min(quantile(psi_rd(~isnan(psi_rd)),[0.025,0.975])) < psi_rd_true) & (psi_rd_true < max(quantile(psi_rd(~isnan(psi_rd)),[0.025,0.975])));

% Calculate width of credible interval
width_psi_rd = abs(max(quantile(psi_rd(~isnan(psi_rd)),[0.025,0.975]))-min(quantile(psi_rd(~isnan(psi_rd)),[0.025,0.975])) );
% end