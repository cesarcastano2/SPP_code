%% plots youngs vs old 
load main_kinem_y
load main_kinem_o

w = 0.1;
% for m=1:length(subjs)
    for i=1:length(conds)
        figure(100)
        subplot(4,2,1)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sl.(conds{i}).avg),mean(kinem_y.sl.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length avg')
% errorbar([1; 2; 3; 4; 5],mean(sl.(conds{i}).avg),mean(sl.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,2)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sl.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sl.(conds{i}).standdev),mean(std(kinem_y.sl.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length std')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,3)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,4)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sf.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sf.(conds{i}).standdev),mean(std(kinem_y.sf.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency std')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,5)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,6)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sw.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sw.(conds{i}).standdev),mean(std(kinem_y.sw.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width std')

for i=1:length(conds)
        figure(100)
        subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,8)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.ws.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i - w,kinem_y.ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.ws.(conds{i}).standdev),mean(std(kinem_y.ws.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')

%old adults 

for i=1:length(conds)
        figure(100)
        subplot(4,2,1)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).avg),mean(kinem_o.sl.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length avg')
% errorbar([1; 2; 3; 4; 5],mean(sl.(conds{i}).avg),mean(sl.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,2)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sl.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).standdev),mean(std(kinem_o.sl.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length std')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,3)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).avg),mean(kinem_o.sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_o.sf.(conds{i}).avg),mean(kinem_o.sf.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,4)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sf.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).standdev),mean(std(kinem_o.kinem_o.sf.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency std')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,5)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).avg),mean(kinem_o.sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_o.sw.(conds{i}).avg),mean(kinem_o.sw.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,6)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.sw.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).standdev),mean(std(kinem_o.sw.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width std')

for i=1:length(conds)
        figure(100)
        subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).avg),mean(kinem_o.ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')
% errorbar([1; 2; 3; 4; 5],mean(kinem_o.ws.(conds{i}).avg),mean(kinem_o.ws.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,8)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i + w,mean(kinem_o.ws.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i + w,kinem_o.ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).standdev),mean(std(kinem_o.ws.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')