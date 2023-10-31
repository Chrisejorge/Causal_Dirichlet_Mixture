function plot_yclusNum_hist(yclusNum)
% Plot the number of occurance of k by histogram

figure
histogram(yclusNum,unique(yclusNum));
annotation('textbox',[.6 .8 .1 .1],'String', [' mean: ' num2str(mean(yclusNum))]);
title('The number of y clusters of data');
xlabel('K');
ylabel('The number of occurance of  k ');
% end