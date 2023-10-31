function plot_psi_rd(psi_rd)
% Plot causal effect

tmp_rd = psi_rd(~isnan(psi_rd));
figure
plot(1:size(tmp_rd,1),tmp_rd,'k');
title('The trace plot of \psi_{rd}');
xlabel('Sample');
ylabel('\psi_{rd}');

% end