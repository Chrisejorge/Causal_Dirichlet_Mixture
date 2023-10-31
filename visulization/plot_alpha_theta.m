function plot_alpha_theta(alpha_theta_draws)
% Plot alpha thetas

figure
plot(1:size(alpha_theta_draws,1),alpha_theta_draws,'k');
title('The trace plot of \alpha_{\theta}');
xlabel('Sample');
ylabel('\alpha_{\theta}');
% end