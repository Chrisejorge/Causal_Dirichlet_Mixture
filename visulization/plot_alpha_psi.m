function plot_alpha_psi(alpha_psi_draws)
% Plot alpha psis

figure
plot(1:size(alpha_psi_draws,1),alpha_psi_draws,'k');
title('The trace plot of \alpha_{\psi}');
xlabel('Sample');
ylabel('\alpha_{\psi}');
% end