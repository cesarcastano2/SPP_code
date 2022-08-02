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
    kinem_y.(kinem_nam{i}).all_var(:,m)=(kinem_y.(kinem_nam{i}).(conds{m}).vardev)';
    end
end

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
          boxchart(z(i,:),kinem_y.sl.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sl.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sl.(conds{i}).standdev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sl.(conds{i}).standdev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sl.(conds{i}).avg),mean(kinem_y.sl.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length std')

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
%         errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).standdev),'k.'), hold on
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
          boxchart(z(i,:),kinem_y.sf.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sf.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sf.(conds{i}).standdev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sf.(conds{i}).standdev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sf.(conds{i}).avg),mean(kinem_y.sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step freq std')

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
%         errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).standdev),'k.'), hold on
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
          boxchart(z(i,:),kinem_y.sw.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.sw.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on

%         plot(i,mean(kinem_y.sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
          plot(x(i)- w,kinem_y.sw.(conds{i}).standdev,'k.','MarkerSize',4), hold on
          plot(x2(i) + w2,kinem_o.sw.(conds{i}).standdev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.sw.(conds{i}).avg),mean(kinem_y.sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width std')

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
%         errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).standdev),'k.'), hold on
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
%           boxplot(z(i,:),kinem_y.ws.(conds{i}).standdev,'Colors','k','Whisker',3,'PlotStyle','compact'), hold on
%           boxplot(z2(i,:),kinem_o.ws.(conds{i}).standdev,'Colors','r','Whisker',3,'PlotStyle','compact'), hold on
          boxchart(z(i,:),kinem_y.ws.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','k'), hold on
          boxchart(z2(i,:),kinem_o.ws.(conds{i}).standdev,'MarkerStyle','none','BoxWidth',0.3,'BoxFaceColor','r'), hold on
%         plot(i,mean(kinem_y.ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
%           plot(x(i)- w,kinem_y.ws.(conds{i}).standdev,'k.','MarkerSize',4), hold on
%           plot(x2(i) + w2,kinem_o.ws.(conds{i}).standdev,'r.','MarkerSize',4), hold on
%         errorbar(i - w,mean(kinem_y.ws.(conds{i}).avg),mean(kinem_y.ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(111)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')
%% detrended old vs young - subplots
load detrend_y
load detrend_o

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];
m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit);
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend);
yng_stack(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_y.(i).speedtrend)];
m = m + 1;
end

figure(112)
subplot(1,4,1)
bar([1,2,3,4,5],yng_stack,'stacked')
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
ylim([0 0.003])
ylabel('variance (m^2)')
title('Young stacked')

m = 1;
for i=conds_s    
detrend_o.(i).det_avg = mean(detrend_o.(i).slminusfit);
detrend_o.(i).speedt_avg = mean(detrend_o.(i).speedtrend);
old_stack(m,:) = [mean(detrend_o.(i).slminusfit), mean(detrend_o.(i).speedtrend)];
m = m + 1;
end

figure(112)
subplot(1,4,2)
bar([1,2,3,4,5],old_stack,'stacked')
ylim([0 0.003])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
ylabel('variance (m^2)')
title('Old stacked')

m = 1;
comb_stack = zeros(10,2);
detrend_y.stack_p=[];
for i=conds_s    
detrend_y.(i).det_avg = mean(detrend_y.(i).slminusfit);
detrend_o.(i).det_avg = mean(detrend_o.(i).slminusfit);
comb_stack_detrend(m,:) = [mean(detrend_y.(i).slminusfit), mean(detrend_o.(i).slminusfit)];
m = m + 1;
end

figure(112)
subplot(1,4,3)
bar([1,2,3,4,5],comb_stack_detrend)
ylim([0 0.002])
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
title('Detrend young vs old')

m = 1;
for i=conds_s    
detrend_y.(i).speedt_avg = mean(detrend_y.(i).speedtrend);
detrend_o.(i).speedt_avg = mean(detrend_o.(i).speedtrend);
comb_stack_speedt(m,:) = [mean(detrend_y.(i).speedtrend), mean(detrend_o.(i).speedtrend)];
m = m + 1;
end

figure(112)
subplot(1,4,4)
bar([1,2,3,4,5],comb_stack_speedt)
ylim([0 0.002]) 
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
title('Speedtrend young vs old')

