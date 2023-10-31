function plot_alpha_omega(alpha_omega_draws)
% Plot alpha omegas

figure
plot(1:size(alpha_omega_draws,1),alpha_omega_draws,'k');
title('The trace plot of \alpha_{\omega}');
xlabel('Sample');
ylabel('\alpha_{\omega}');
% end