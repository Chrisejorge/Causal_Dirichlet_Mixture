function plot_psi_rr(psi_rr)
% Plot causal relative risk

tmp_rr = psi_rr(~isnan(psi_rr));
figure
plot(1:size(tmp_rr,1),tmp_rr,'k');
title('The trace plot of \psi_{rr}');
xlabel('Sample');
ylabel('\psi_{rr}');

% end