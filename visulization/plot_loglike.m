function plot_loglike(loglikeregData)
% Plot log likelihood of data

figure
plot(1:1:size(loglikeregData,1),loglikeregData,'k');
title('The log likelihood of data');
xlabel('Sample');
ylabel('Log likelihood');
legend('Log likelihood');
% end