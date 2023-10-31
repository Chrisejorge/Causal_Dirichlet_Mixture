function plot_yclusNum(yclusNum)
% Plot the number of occurance of k

figure
plot(1:size(yclusNum,1),yclusNum,'k');
title('The number of y clusters of data');
xlabel('Sample');
ylabel('K');
% end