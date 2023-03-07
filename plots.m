%% plots youngs vs old 
load main_kinem_y
load main_kinem_o
load detrend_y
load detrend_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};

%% old vs young kinematic with individual subjects
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).avg),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).avg),mean(kinem_o.sl.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).standdev),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).standdev),mean(std(kinem_o.sl.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).avg),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).avg),mean(kinem_o.sf.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).standdev),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).standdev),mean(std(kinem_o.sf.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).avg),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).avg),mean(kinem_o.sw.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).standdev),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).standdev),mean(std(kinem_o.sw.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).avg),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).avg),mean(kinem_o.ws.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).standdev),'r.','MarkerSize',20), hold on
        plot(i + w,kinem_o.ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).standdev),mean(std(kinem_o.ws.(conds{i}).standdev)),'r.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')
%% plots youngs vs old - no individual subjects 
load main_kinem_y
load main_kinem_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
w = 0.1;
% for m=1:length(subjs)
    for i=1:length(conds)
        figure(100)
        subplot(4,2,1)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i - w,mean(kinem_y.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
%         plot(i - w,kinem_y.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
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
%         plot(i - w,kinem_y.ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).avg),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).avg),mean(kinem_o.sl.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).standdev),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).standdev),mean(std(kinem_o.sl.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).avg),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).avg),mean(kinem_o.sf.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).standdev),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).standdev),mean(std(kinem_o.sf.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).avg),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).avg),mean(kinem_o.sw.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).standdev),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).standdev),mean(std(kinem_o.sw.(conds{i}).standdev)),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).avg),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).avg),mean(kinem_o.ws.(conds{i}).standdev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).standdev),'r.','MarkerSize',20), hold on
%         plot(i + w,kinem_o.ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).standdev),mean(std(kinem_o.ws.(conds{i}).standdev)),'r.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')

%% detrended old vs young
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];
m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit);
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend);
comb_stack(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_y.(i).speedtrend)];
m = m + 2;
end
m = 2;
for i=conds_s    
detrend_o.(i).det_avg = mean(detrend_o.(i).slminusfit);
detrend_o.(i).speedt_avg = mean(detrend_o.(i).speedtrend);
comb_stack(m,:) = [mean(detrend_o.(i).slminusfit), mean(detrend_o.(i).speedtrend)];
m = m + 2;
end

figure
bar([0.9,1.1,1.9,2.1,2.9,3.1,3.9,4.1,4.9,5.1],comb_stack,'stacked')
% legend('detrended','speedtrend')
% set(gca, 'XTick', [0.1,1.1,1.9,2.1,2.9,3.1,3.9,4.1,4.9,5.1,5.9,6.1,6.9,7.1,7.9,8.1,8.9,9.1,9.9,10.1],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
ylabel('variance (m^2)')
%% rearange the data 
kinem_nam = {'sl' 'sf' 'sw' 'ws'};
conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
for i=1:length(kinem_nam)
    for m=1:length(conds)
    kinem_o.(kinem_nam{i}).all_ws(:,m)=(kinem.(kinem_nam{i}).(conds{m}).avg)';
    end
end
% rearange the data 

kinem_nam = {'sl' 'sf' 'sw' 'ws'};
conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
% for i=1:length(kinem_nam)
    for m=1:length(conds)
    detrend_y.all_det(:,m)=(detrend_y.(conds{m}).slminusfit_vd);
    detrend_y.all_speedt(:,m)=(detrend_y.(conds{m}).speedtrend_vd);
    detrend_o.all_det(:,m)=(detrend_o.(conds{m}).slminusfit_vd);
    detrend_o.all_speedt(:,m)=(detrend_o.(conds{m}).speedtrend_vd);
%     detrend_y.all_totalvar(:,m)=(detrend_y.(conds{m}).totalvar);
%     detrend_o.all_totalvar(:,m)=(detrend_o.(conds{m}).totalvar); 
    end
% end totalvar
%% plots youngs vs old - no individual subjects - box plots 
load main_kinem_y
load main_kinem_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
w = [0.05,0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];
w2 = [0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];

x = [0.8;1.8;2.8;3.8;4.8];
x2 = [1.2;2.2;3.2;4.2;5.2];
z = repmat(x,1,10);
z2 = repmat(x2,1,9);
% for m=1:length(subjs)
    for i=1:length(conds)
        figure(110)
        subplot(3,2,1)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sl.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sl.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sl.(conds{i}).avg,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sl.(conds{i}).avg,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sl.(conds{i}).avg),mean(kinem_y.sl.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length avg')

    for i=1:length(conds)
        figure(110)
        subplot(3,2,2)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sl.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sl.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sl.(conds{i}).vardev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sl.(conds{i}).vardev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sl.(conds{i}).avg),mean(kinem_y.sl.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length var')

for i=1:length(conds)
        figure(110)
        subplot(3,2,3)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sf.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sf.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sf.(conds{i}).avg,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sf.(conds{i}).avg,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step freq avg')

    for i=1:length(conds)
        figure(110)
        subplot(3,2,4)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sf.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sf.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sf.(conds{i}).vardev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sf.(conds{i}).vardev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step freq var')

for i=1:length(conds)
        figure(110)
        subplot(3,2,5)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sw.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sw.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sw.(conds{i}).avg,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sw.(conds{i}).avg,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')

    for i=1:length(conds)
        figure(110)
        subplot(3,2,6)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.sw.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sw.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sw.(conds{i}).vardev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sw.(conds{i}).vardev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width var')

for i=1:length(conds)
        figure(111)
        subplot(1,2,1)
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxchart(z(i,:),kinem_y.ws.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.ws.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.ws.(conds{i}).avg,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.ws.(conds{i}).avg,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).vardev),'k.'), hold on
    end
% end
figure(111)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')

    for i=1:length(conds)
        figure(111)
        subplot(1,2,2)
%         xlim([0 6])
%         ylim([0.6 0.85])
%           boxplot(z(i,:),kinem_y.ws.(conds{i}).vardev,'Colors','k','Whisker',3,'PlotStyle','compact'), hold on
%           boxplot(z2(i,:),kinem_o.ws.(conds{i}).vardev,'Colors','r','Whisker',3,'PlotStyle','compact'), hold on
          boxchart(z(i,:),kinem_y.ws.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.ws.(conds{i}).vardev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on
%         plot(i,mean(kinem_y.ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.ws.(conds{i}).vardev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.ws.(conds{i}).vardev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(111)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed var')
%% detrended SL old vs young - subplots
% load detrend_y
% load detrend_o
conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
detrend_y.det_mean=[]; detrend_y.det_std=[]; detrend_y.speedt_mean=[]; detrend_y.speedt_std=[];
for i=conds_s   
detrend_y.det_mean(m) = mean(detrend_y.(i).slminusfit);
detrend_y.speedt_mean(m) = mean(detrend_y.(i).speedtrend);
detrend_y.det_std(m) = std(detrend_y.(i).slminusfit);
detrend_y.speedt_std(m) = std(detrend_y.(i).speedtrend);
yng_stack(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_y.(i).speedtrend)];
m = m + 1;
end
figure(112)
subplot(1,3,1)
bar([1,2,3,4,5],yng_stack,'stacked'); hold on
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_y.det_mean,detrend_y.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrend_y.speedt_mean + detrend_y.det_mean),detrend_y.speedt_std,'.');
ylim([0 0.006])
ylabel('variance (m^2)')
title('Young SL stacked')

m = 1;
detrend_o.stack_p=[];
detrend_o.det_mean=[]; detrend_o.det_std=[]; detrend_o.speedt_mean=[]; detrend_o.speedt_std=[];
for i=conds_s   
detrend_o.det_mean(m) = mean(detrend_o.(i).slminusfit);
detrend_o.speedt_mean(m) = mean(detrend_o.(i).speedtrend);
detrend_o.det_std(m) = std(detrend_o.(i).slminusfit);
detrend_o.speedt_std(m) = std(detrend_o.(i).speedtrend);
old_stack(m,:) = [mean(detrend_o.(i).slminusfit), mean(detrend_o.(i).speedtrend)];
m = m + 1;
end
figure(112)
subplot(1,3,2)
bar([1,2,3,4,5],old_stack,'stacked'); hold on;
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_o.det_mean,detrend_o.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrend_o.speedt_mean + detrend_o.det_mean),detrend_o.speedt_std,'.');
ylabel('variance (m^2)')
title('Old SL stacked')

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit);
detrend_o.(i).det_avg = mean(detrend_o.(i).slminusfit);
comb_stack_detrend(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_o.(i).slminusfit)];
m = m + 1;
end
x_det_yng = [0.85,1.85,2.85,3.85,4.85];
x_det_old = [1.15,2.15,3.15,4.15,5.15];
figure(112)
subplot(2,3,3)
bar([1,2,3,4,5],comb_stack_detrend); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrend_y.det_mean,detrend_y.det_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrend_o.det_mean,detrend_o.det_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).slminusfit,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrend_o.(i).slminusfit,'.','MarkerSize',4), hold on
m = m + 1;
end
title('Detrend SL young vs old')

m = 1;
for i=conds_s    
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend);
detrend_o.(i).speedt_avg = mean(detrend_o.(i).speedtrend);
comb_stack_speedt(m,:) = [mean(detrend_y.(i).speedtrend), mean(detrend_o.(i).speedtrend)];
m = m + 1;
end

figure(112)
subplot(2,3,6)
bar([1,2,3,4,5],comb_stack_speedt); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrend_y.speedt_mean,detrend_y.speedt_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrend_o.speedt_mean,detrend_o.speedt_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).speedtrend,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrend_o.(i).speedtrend,'.','MarkerSize',4), hold on
m = m + 1;
end
title('Speedtrend SL young vs old')
%% detrended step width old vs young - subplots
load detrend_w_y
load detrend_w_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];

m = 1;
comb_stack = zeros(10,2);
detrended_w_y.stack_p=[];
detrended_w_y.det_mean=[]; detrended_w_y.det_std=[]; detrended_w_y.speedt_mean=[]; detrended_w_y.speedt_std=[];
for i=conds_s   
detrended_w_y.det_mean(m) = mean(detrended_w_y.(i).slminusfit);
detrended_w_y.speedt_mean(m) = mean(detrended_w_y.(i).speedtrend);
detrended_w_y.det_std(m) = std(detrended_w_y.(i).slminusfit);
detrended_w_y.speedt_std(m) = std(detrended_w_y.(i).speedtrend);
yng_stack(m,:) = [mean(detrended_w_y.(i).slminusfit), mean(detrended_w_y.(i).speedtrend)];
m = m + 1;
end
figure(113)
subplot(1,3,1)
bar([1,2,3,4,5],yng_stack,'stacked'); hold on
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrended_w_y.det_mean,detrended_w_y.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrended_w_y.speedt_mean + detrended_w_y.det_mean),detrended_w_y.speedt_std,'.');
ylim([0 0.002])
ylabel('variance (m^2)')
title('Young SW stacked')

m = 1;
detrended_w_o.stack_p=[];
detrended_w_o.det_mean=[]; detrended_w_o.det_std=[]; detrended_w_o.speedt_mean=[]; detrended_w_o.speedt_std=[];
for i=conds_s   
detrended_w_o.det_mean(m) = mean(detrended_w_o.(i).slminusfit);
detrended_w_o.speedt_mean(m) = mean(detrended_w_o.(i).speedtrend);
detrended_w_o.det_std(m) = std(detrended_w_o.(i).slminusfit);
detrended_w_o.speedt_std(m) = std(detrended_w_o.(i).speedtrend);
old_stack(m,:) = [mean(detrended_w_o.(i).slminusfit), mean(detrended_w_o.(i).speedtrend)];
m = m + 1;
end
figure(113)
subplot(1,3,2)
bar([1,2,3,4,5],old_stack,'stacked'); hold on;
ylim([0 0.002])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrended_w_o.det_mean,detrended_w_o.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrended_w_o.speedt_mean + detrended_w_o.det_mean),detrended_w_o.speedt_std,'.');
ylabel('variance (m^2)')
title('Old SW stacked')

m = 1;
comb_stack = zeros(10,2);
detrended_w_y.stack_p=[];
for i=conds_s    
detrended_w_y.(i).det_avg = mean(detrended_w_y.(i).slminusfit);
detrended_w_o.(i).det_avg = mean(detrended_w_o.(i).slminusfit);
comb_stack_detrend(m,:) = [mean(detrended_w_y.(i).slminusfit), mean(detrended_w_o.(i).slminusfit)];
m = m + 1;
end
x_det_yng = [0.85,1.85,2.85,3.85,4.85];
x_det_old = [1.15,2.15,3.15,4.15,5.15];
figure(113)
subplot(2,3,3)
bar([1,2,3,4,5],comb_stack_detrend); hold on
ylim([0 0.002])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrended_w_y.det_mean,detrended_w_y.det_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrended_w_o.det_mean,detrended_w_o.det_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrended_w_y.(i).slminusfit,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrended_w_o.(i).slminusfit,'.','MarkerSize',4), hold on
m = m + 1;
end
title('Detrend SW young vs old')

m = 1;
for i=conds_s    
detrended_w_y.(i).speedt_avg = mean(detrended_w_y.(i).speedtrend);
detrended_w_o.(i).speedt_avg = mean(detrended_w_o.(i).speedtrend);
comb_stack_speedt(m,:) = [mean(detrended_w_y.(i).speedtrend), mean(detrended_w_o.(i).speedtrend)];
m = m + 1;
end

figure(113)
subplot(2,3,6)
bar([1,2,3,4,5],comb_stack_speedt); hold on
ylim([0 0.002])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrended_w_y.speedt_mean,detrended_w_y.speedt_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrended_w_o.speedt_mean,detrended_w_o.speedt_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrended_w_y.(i).speedtrend,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrended_w_o.(i).speedtrend,'.','MarkerSize',4), hold on
m = m + 1;
end
title('Speedtrend SW young vs old')

%% detrended SL old vs young - subplots
load detrend_y
load detrend_o
conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
detrend_y.det_mean=[]; detrend_y.det_std=[]; detrend_y.speedt_mean=[]; detrend_y.speedt_std=[];
for i=conds_s   
detrend_y.det_mean(m) = mean(detrend_y.(i).slminusfit);
detrend_y.speedt_mean(m) = mean(detrend_y.(i).speedtrend);
detrend_y.det_std(m) = std(detrend_y.(i).slminusfit);
detrend_y.speedt_std(m) = std(detrend_y.(i).speedtrend);
yng_stack(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_y.(i).speedtrend)];
m = m + 1;
end
figure(112)
subplot(1,2,1)
bar([1,2,3,4,5],yng_stack,'stacked'); hold on
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_y.det_mean,detrend_y.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrend_y.speedt_mean + detrend_y.det_mean),detrend_y.speedt_std,'.');
ylim([0 0.006])
ylabel('variance (m^2)')
title('Young SL stacked')

%% abstract
load main_kinem_y
subjs = {'SPP2' 'SPP3' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11'};

color = jet(length(subjs)); 
for i=1:length(conds)
        figure(199)
%         subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
         for m=1:4
            plot(i+(0.02*m),ws.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
         end
        for m=5:8
            
            plot(i+(-0.02*m),ws.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(ws.(conds{i}).avg),mean(ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(199)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')

for i=1:length(conds)
        figure(198)
%         subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
         for m=1:4
            plot(i+(0.02*m),sw.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
         end
        for m=5:8
            
            plot(i+(-0.02*m),sw.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sw.(conds{i}).avg),mean(sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(198)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')

for i=1:length(conds)
        figure(197)
%         subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
         for m=1:4
            plot(i+(0.02*m),sf.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
         end
        for m=5:8
            
            plot(i+(-0.02*m),sf.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sf.(conds{i}).avg),mean(sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(197)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency avg')

for i=1:length(subjs)
    speed_spss(i,:) = ws.([subjs{i}]).avg;
end 
for i=1:length(subjs)
    frequency_spss(i,:) = sf.([subjs{i}]).avg;
end 
for i=1:length(subjs)
    width_spss(i,:) = sw.([subjs{i}]).avg;
end 


m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit);
comb_stack_detrend(m,:) = [mean(detrend_y.(i).det_avg)];
m = m + 1;
end
x_det_yng = [1,2,3,4,5];
figure(112)
subplot(2,2,2)
bar([1,2,3,4,5],comb_stack_detrend); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_y.det_mean,detrend_y.det_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).slminusfit,'.','MarkerSize',10), hold on
m = m + 1;
end
title('Detrend SL young vs old')

m = 1;
for i=conds_s    
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend);
comb_stack_speedt(m,:) = [mean(detrend_y.(i).speedtrend)];
m = m + 1;
end

figure(112)
subplot(2,2,4)
bar([1,2,3,4,5],comb_stack_speedt); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_y.speedt_mean,detrend_y.speedt_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).speedtrend,'.','MarkerSize',10), hold on
m = m + 1;
end
title('Speedtrend SL young vs old')


%% NEW BOX PLOTS
load main_kinem_y
load main_kinem_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
w = [0.05,0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];
w2 = [0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];

x = [0.8;1.8;2.8;3.8;4.8];
x2 = [1.2;2.2;3.2;4.2;5.2];
z = repmat(x,1,10);
z2 = repmat(x2,1,9);

boxp.one = boxplot(boxplots.steplengthvar_speedtrend.plot_values,'Whisker',3);
for i=1:length(conds)
       
%         xlim([0 6])
%         ylim([0.6 0.85])
          boxpvy.sl(:,i) = kinem_y.sl.(conds{i}).avg' 
          boxpvo.sl(:,i) = kinem_o.sl.(conds{i}).avg'
          boxpvy.sf(:,i) = kinem_y.sf.(conds{i}).avg' 
          boxpvo.sf(:,i) = kinem_o.sf.(conds{i}).avg'
          boxpvy.sw(:,i) = kinem_y.sw.(conds{i}).avg' 
          boxpvo.sw(:,i) = kinem_o.sw.(conds{i}).avg'
          boxpvy.ws(:,i) = kinem_y.ws.(conds{i}).avg' 
          boxpvo.ws(:,i) = kinem_o.ws.(conds{i}).avg'
          boxpvy.slvar(:,i) = kinem_y.sl.(conds{i}).vardev' 
          boxpvo.slvar(:,i) = kinem_o.sl.(conds{i}).vardev'
          boxpvy.sfvar(:,i) = kinem_y.sf.(conds{i}).vardev' 
          boxpvo.sfvar(:,i) = kinem_o.sf.(conds{i}).vardev'
          boxpvy.swvar(:,i) = kinem_y.sw.(conds{i}).vardev' 
          boxpvo.swvar(:,i) = kinem_o.sw.(conds{i}).vardev'
          boxpvy.wsvar(:,i) = kinem_y.ws.(conds{i}).vardev' 
          boxpvo.wsvar(:,i) = kinem_o.ws.(conds{i}).vardev'
%           boxchart(z2(i,:),kinem_o.sf.(conds{i}).avg,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
%           plot(x(i)- w,kinem_y.sf.(conds{i}).avg,'k.','MarkerSize',4), hold on
%           plot(x2(i) + w2,kinem_o.sf.(conds{i}).avg,'r.','MarkerSize',4), hold on
% %         errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).vardev),'k.'), hold on
end

figure(309)
boxplot(boxpvy.sw,'Whisker',3); hold on
for i=1:length(conds)
plot(i,boxpvy.sw(:,i),'.','MarkerSize',10), hold on
end 
figure(310)
boxplot(boxpvo.sl,'Whisker',3); hold on
for i=1:length(conds)
plot(i,boxpvo.sl(:,i),'.','MarkerSize',10), hold on
end 
figure(311)
boxplot(boxpvo.slvar,'Whisker',5); hold on
for i=1:length(conds)
plot(i,boxpvo.slvar(:,i),'.','MarkerSize',10), hold on
end 
figure(312)
boxplot(boxpvo.sfvar,'Whisker',10); hold on
for i=1:length(conds)
plot(i,boxpvo.sfvar(:,i),'.','MarkerSize',10), hold on
end 

title('SL AVG young')
figure(310)
boxplot(boxpvy.sw,'Whisker',3);
title('SW AVG young')
figure(311)
boxplot(boxpvo.sw,'Whisker',3);
title('SW AVG old')

%% old vs young kinematic with individual subjects
% x = [0.8;1.8;2.8;3.8;4.8];
% x2 = [1.2;2.2;3.2;4.2;5.2];
w = 0.1;
wy = [0.05,0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];
wo = [0.054,0.058,0.062,0.066,0.07,0.074,0.076,0.08,0.084];
wy = wy*(4);
wo = wo*(4);
% for m=1:length(subjs)
    for i=1:length(conds)
        figure(100)
        subplot(4,2,1)
        xlim([0 6])
        ylim([0.583785901645933 0.832296864500448])
        plot(i - w,mean(kinem_y.sl.(conds{i}).avg),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sl.(conds{i}).avg),std(kinem_y.sl.(conds{i}).avg),'b.'), hold on
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
        ylim([0 0.006])
        plot(i - w,mean(kinem_y.sl.(conds{i}).vardev),'b.','MarkerSize',10), hold on
%         plot(i - wy,kinem_y.sl.(conds{i}).vardev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sl.(conds{i}).vardev), std(kinem_y.sl.(conds{i}).vardev),'b.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length var')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,3)
        xlim([0 6])
        ylim([1.739936389098314 2.5])
        plot(i - w,mean(kinem_y.sf.(conds{i}).avg),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),std(kinem_y.sf.(conds{i}).avg),'b.'), hold on
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
        ylim([0 0.0582287746936])
        plot(i - w,mean(kinem_y.sf.(conds{i}).vardev),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.sf.(conds{i}).vardev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sf.(conds{i}).vardev),std(kinem_y.sf.(conds{i}).vardev),'b.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency var')

    for i=1:length(conds)
        figure(100)
        subplot(4,2,5)
        xlim([0 6])
        ylim([0.069994841906375 0.2])
        plot(i - w,mean(kinem_y.sw.(conds{i}).avg),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),std(kinem_y.sw.(conds{i}).avg),'b.'), hold on
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
        ylim([0 0.002139613829304])
        plot(i - w,mean(kinem_y.sw.(conds{i}).vardev),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.sw.(conds{i}).vardev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.sw.(conds{i}).vardev),std(kinem_y.sw.(conds{i}).vardev),'b.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width var')

for i=1:length(conds)
        figure(100)
        subplot(4,2,7)
        xlim([0 6])
        ylim([1.1 1.9])
        plot(i - w,mean(kinem_y.ws.(conds{i}).avg),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),std(kinem_y.ws.(conds{i}).avg),'b.'), hold on
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
        ylim([0 0.04])
        plot(i - w,mean(kinem_y.ws.(conds{i}).vardev),'b.','MarkerSize',10), hold on
        plot(i - wy,kinem_y.ws.(conds{i}).vardev,'b.','MarkerSize',4), hold on
        errorbar(i - w,mean(kinem_y.ws.(conds{i}).vardev),std(kinem_y.ws.(conds{i}).vardev),'b.'), hold on
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).avg),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.sl.(conds{i}).avg,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).avg),std(kinem_o.sl.(conds{i}).avg),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sl.(conds{i}).vardev),'r.','MarkerSize',10), hold on
%         plot(i + wo,kinem_o.sl.(conds{i}).vardev,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sl.(conds{i}).vardev),std(kinem_o.sl.(conds{i}).vardev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).avg),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.sf.(conds{i}).avg,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).avg),std(kinem_o.sf.(conds{i}).avg),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sf.(conds{i}).vardev),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.sf.(conds{i}).vardev,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sf.(conds{i}).vardev),std(kinem_o.sf.(conds{i}).vardev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).avg),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.sw.(conds{i}).avg,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).avg),std(kinem_o.sw.(conds{i}).avg),'r.'), hold on
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
        plot(i + w,mean(kinem_o.sw.(conds{i}).vardev),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.sw.(conds{i}).vardev,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.sw.(conds{i}).vardev),std(kinem_o.sw.(conds{i}).vardev),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).avg),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.ws.(conds{i}).avg,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).avg),std(kinem_o.ws.(conds{i}).avg),'r.'), hold on
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
        plot(i + w,mean(kinem_o.ws.(conds{i}).vardev),'r.','MarkerSize',10), hold on
        plot(i + wo,kinem_o.ws.(conds{i}).vardev,'r.','MarkerSize',4), hold on
        errorbar(i + w,mean(kinem_o.ws.(conds{i}).vardev),std(kinem_o.ws.(conds{i}).vardev),'r.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')

%% detrended total var 
load main_kinem_y
load main_kinem_o
load detrend_y
load detrend_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};

for i=1:length(conds)
detrend_y.(conds{i}).vardiff = (kinem_y.sl.(conds{i}).vardev')./detrend_y.(conds{i}).totalvar;
detrend_y.(conds{i}).slminusfit_vd = detrend_y.(conds{i}).vardiff .* detrend_y.(conds{i}).slminusfit;
detrend_y.(conds{i}).speedtrend_vd = detrend_y.(conds{i}).vardiff .* detrend_y.(conds{i}).speedtrend;
detrend_y.(conds{i}).vardiff_sum = detrend_y.(conds{i}).slminusfit_vd + detrend_y.(conds{i}).speedtrend_vd;
detrend_y.(conds{i}).totaldiff = detrend_y.(conds{i}).vardiff_sum - (kinem_y.sl.(conds{i}).vardev');
detrend_o.(conds{i}).vardiff = (kinem_o.sl.(conds{i}).vardev')./detrend_o.(conds{i}).totalvar;
detrend_o.(conds{i}).slminusfit_vd = detrend_o.(conds{i}).vardiff .* detrend_o.(conds{i}).slminusfit;
detrend_o.(conds{i}).speedtrend_vd = detrend_o.(conds{i}).vardiff .* detrend_o.(conds{i}).speedtrend;
detrend_o.(conds{i}).vardiff_sum = detrend_o.(conds{i}).slminusfit_vd + detrend_o.(conds{i}).speedtrend_vd;
detrend_o.(conds{i}).totaldiff = detrend_o.(conds{i}).vardiff_sum - (kinem_o.sl.(conds{i}).vardev');
end 

%% detrended SL old vs young - subplots
% load detrend_y
% load detrend_o
conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
detrend_y.det_mean=[]; detrend_y.det_std=[]; detrend_y.speedt_mean=[]; detrend_y.speedt_std=[];
for i=conds_s   
detrend_y.det_mean(m) = mean(detrend_y.(i).slminusfit_vd);
detrend_y.speedt_mean(m) = mean(detrend_y.(i).speedtrend_vd);
detrend_y.det_std(m) = std(detrend_y.(i).slminusfit_vd);
detrend_y.speedt_std(m) = std(detrend_y.(i).speedtrend_vd);
yng_stack(m,:) = [mean(detrend_y.(i).slminusfit_vd), mean(detrend_y.(i).speedtrend_vd)];
m = m + 1;
end
figure(112)
subplot(1,3,1)
bar([1,2,3,4,5],yng_stack,'stacked'); hold on
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_y.det_mean,(detrend_y.det_std - detrend_y.det_std),detrend_y.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrend_y.speedt_mean + detrend_y.det_mean),(detrend_y.speedt_std - detrend_y.speedt_std),detrend_y.speedt_std,'.');
ylim([0 0.006])
ylabel('variance (m^2)')
title('Young SL stacked')

m = 1;
detrend_o.stack_p=[];
detrend_o.det_mean=[]; detrend_o.det_std=[]; detrend_o.speedt_mean=[]; detrend_o.speedt_std=[];
for i=conds_s   
detrend_o.det_mean(m) = mean(detrend_o.(i).slminusfit_vd);
detrend_o.speedt_mean(m) = mean(detrend_o.(i).speedtrend_vd);
detrend_o.det_std(m) = std(detrend_o.(i).slminusfit_vd);
detrend_o.speedt_std(m) = std(detrend_o.(i).speedtrend_vd);
old_stack(m,:) = [mean(detrend_o.(i).slminusfit_vd), mean(detrend_o.(i).speedtrend_vd)];
m = m + 1;
end
figure(112)
subplot(1,3,2)
bar([1,2,3,4,5],old_stack,'stacked'); hold on;
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([1,2,3,4,5],detrend_o.det_mean,(detrend_o.det_std - detrend_o.det_std) ,detrend_o.det_std,'.'); hold on;
errorbar([1,2,3,4,5],(detrend_o.speedt_mean + detrend_o.det_mean),(detrend_o.speedt_std - detrend_o.speedt_std) ,detrend_o.speedt_std,'.');
ylabel('variance (m^2)')
title('Old SL stacked')

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit_vd);
detrend_o.(i).det_avg = mean(detrend_o.(i).slminusfit_vd);
comb_stack_detrend(m,:) = [mean(detrend_y.(i).slminusfit_vd), mean(detrend_o.(i).slminusfit_vd)];
m = m + 1;
end
x_det_yng = [0.85,1.85,2.85,3.85,4.85];
x_det_old = [1.15,2.15,3.15,4.15,5.15];
figure(112)
subplot(2,3,3)
bar([1,2,3,4,5],comb_stack_detrend); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrend_y.det_mean,(detrend_y.det_std - detrend_y.det_std) ,detrend_y.det_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrend_o.det_mean,(detrend_o.det_std - detrend_o.det_std) ,detrend_o.det_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).slminusfit_vd,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrend_o.(i).slminusfit_vd,'.','MarkerSize',4), hold on
m = m + 1;
end
title('Detrend SL young vs old')

m = 1;
for i=conds_s    
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend_vd);
detrend_o.(i).speedt_avg = mean(detrend_o.(i).speedtrend_vd);
comb_stack_speedt(m,:) = [mean(detrend_y.(i).speedtrend_vd), mean(detrend_o.(i).speedtrend_vd)];
m = m + 1;
end

figure(112)
subplot(2,3,6)
bar([1,2,3,4,5],comb_stack_speedt); hold on
ylim([0 0.006])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
errorbar([0.85,1.85,2.85,3.85,4.85],detrend_y.speedt_mean,(detrend_y.speedt_std - detrend_y.speedt_std) ,detrend_y.speedt_std,'.'); hold on;
errorbar([1.15,2.15,3.15,4.15,5.15],detrend_o.speedt_mean,(detrend_o.speedt_std - detrend_o.speedt_std) ,detrend_o.speedt_std,'.'); hold on;
m = 1;
for i=conds_s
plot(x_det_yng(m),detrend_y.(i).speedtrend_vd,'.','MarkerSize',4), hold on
plot(x_det_old(m),detrend_o.(i).speedtrend_vd,'.','MarkerSize',4), hold on
m = m + 1;
end
title('speedtrend_vd SL young vs old')