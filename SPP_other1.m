%% NEW CODE 
close all
clearvars
clc

% subjs = {'SPP2' 'SPP3' 'SPP4' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11' 'SPP12'};
% subjs = {'SPP13' 'SPP14' 'SPP15' 'SPP16' 'SPP17' 'SPP19' 'SPP21' 'SPP22' 'SPP24'};
% subjs = {'SPP2' 'SPP4' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11' 'SPP12'};
subjs = {'SPP13' 'SPP14' 'SPP15' 'SPP17' 'SPP19' 'SPP21' 'SPP22' 'SPP24'};

% subjs = {'SPP3'};

conds_f = {'0' '1' '2' '3' '4'};
% conds_f = {'4'};

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
% conds = {'diff_fm'};


fs = 240; %hz, resamp tm to df rate
proj = 'F:\SPP\subjects\';
% proj = 'Z:\SPP\subjects\';
% proj = 'X:\SPP\subjects'; %NAS COM_allputer


for s = 1:length(subjs)
    for c = 1:length(conds_f)
       dtm=[]; forces_df=[]; GEgood=[]; GE=[]; sl_temp=[]; steplength_speed=[]; steplength_time=[]; sw_temp=[];
        
        dflow_file = [proj subjs{s} '\perturb_' conds_f{c} '0001.txt'];
        dflow_treadmill_file = [proj subjs{s} '\perturb_' conds_f{c} '_treadmill0001.txt'];
       
% dflow_file = ['F:\SPP\Helen\S20_diff_f0001.txt'];

          if ~exist('dflow_file','var')
              [FileName,PathName,FilterIndex]  = uigetfile('*.txt');
              dflow_file = [PathName FileName];
          end
          
          if ~exist('dflow_treadmill_file','var')
            [FileName,PathName,FilterIndex]  = uigetfile('*treadmill0001.txt');
            dflow_treadmill_file = [PathName FileName];
          end
          
          [Frame_df, Time_df, markers_df, forces_df, startidx, stopidx, Total] = import_dflow(dflow_file,subjs(s),conds_f(c));
         
          datatreadmill_all.(conds{c}) = importTreadmillFile_SPP2(dflow_treadmill_file);
%           
%           full_data_treadmill = dtm; 


          for i=1:length(markers_df.labels)
           markers_df.(markers_df.labels{i}) = markers_df.(markers_df.labels{i});
          end
          
%            treadmill_file = [subjs{s} '_' conds{c} 't'];
           
%            tf.(conds{c}) = eval(treadmill_file);
            
        tf.speed.(conds{c})= datatreadmill_all.(conds{c}).SpeedActual1;
        tf.sway.(conds{c})= datatreadmill_all.(conds{c}).SwayActual;
        tf.stride.(conds{c})= datatreadmill_all.(conds{c}).stride;
        tf.time.(conds{c})= datatreadmill_all.(conds{c}).Time;
        tf.time_mark.(conds{c})= datatreadmill_all.(conds{c}).t;
        
  
%% adjust resfresh rates 
       dtm.Time = interpft(tf.time.(conds{c})(:),length(Frame_df));
       dtm.LeftBeltSpeed = interpft(tf.speed.(conds{c})(:),length(Frame_df));
       dtm.RightBeltSpeed  = interpft(tf.speed.(conds{c})(:),length(Frame_df));
       dtm.Time_mark = interpft(tf.time_mark.(conds{c})(:),length(Frame_df));
       dtm.stride = interpft(tf.stride.(conds{c})(:),length(Frame_df));
%        dtm.sway2 = interp1(tf.sway.(conds{c})(:),1:length(Frame_df));
       dtm.sway = interpft(tf.sway.(conds{c})(:),length(Frame_df));


       dtm.sway_v1 = (dtm.sway(startidx:stopidx));       
       dtm.stride_v1 = dtm.stride(startidx:stopidx);
       dtm.Time_mark_v1 = dtm.Time_mark(startidx:stopidx);
       dtm.Time_v1 = dtm.Time(startidx:stopidx);
       dtm.LeftBeltSpeed_v1 =  dtm.LeftBeltSpeed(startidx:stopidx);
       dtm.RightBeltSpeed_v1  = dtm.RightBeltSpeed(startidx:stopidx);
       
       for i=1:length(markers_df.labels)
           markers_df_idx.(markers_df.labels{i}) = markers_df.(markers_df.labels{i})(startidx:stopidx,:);
       end
       
       forces_df.FP1Cop_v1 = forces_df.FP1Cop(startidx:stopidx,:);
       forces_df.FP1For_v1 = forces_df.FP1For(startidx:stopidx,:);
       forces_df.FP1Mom_v1 = forces_df.FP1Mom(startidx:stopidx,:);
       forces_df.FP2Cop_v1 = forces_df.FP2Cop(startidx:stopidx,:);
       forces_df.FP2For_v1 = forces_df.FP2For(startidx:stopidx,:);
       forces_df.FP2Mom_v1 = forces_df.FP2Mom(startidx:stopidx,:);
        
       Frame_df_v1 = Frame_df(startidx:stopidx,:);
       Time_df_v1 = Time_df(startidx:stopidx,:);
       

%% start times 
        Dflow_size = length(Time_df);
        
if c==1
    istart = find(dtm.Time_mark_v1 < 30,1,'last');
    Stop_Time = find(dtm.Time_mark_v1 < (round(dtm.Time_mark_v1(end))-2),1,'last');
else
    istart = find(dtm.Time_mark_v1 < 30,1,'last');
    Stop_Time = find(dtm.sway_v1 < -0.005,1,'first');
end
       
        start_values = istart;       
        start_value_pref = start_values;
        Start_Time = start_value_pref;
      

       dtm.Time_mark_s = dtm.Time_mark_v1(Start_Time:Stop_Time,:);
       dtm.Time_s = dtm.Time_v1(Start_Time:Stop_Time,:);
       dtm.LeftBeltSpeed_s =  dtm.LeftBeltSpeed_v1(Start_Time:Stop_Time,:);
       dtm.RightBeltSpeed_s  = dtm.RightBeltSpeed_v1(Start_Time:Stop_Time,:);
       dtm.sway_s = dtm.sway_v1(Start_Time:Stop_Time,:);
         
       
        for i=1:length(markers_df.labels)
            markers_df_s.(markers_df.labels{i}) = markers_df_idx.(markers_df.labels{i})(Start_Time:Stop_Time,:);
        end
        markers_df_s.labels = markers_df.labels;
        
        forces_df.forces_df_s.FP1Cop = forces_df.FP1Cop_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.FP1For = forces_df.FP1For_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.FP1Mom = forces_df.FP1Mom_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.FP2Cop = forces_df.FP2Cop_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.FP2For = forces_df.FP2For_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.FP2Mom = forces_df.FP2Mom_v1(Start_Time:Stop_Time,:);
        forces_df.forces_df_s.labels = forces_df.labels;
        
        Frame_df_s = Frame_df_v1(Start_Time:Stop_Time,:);
        Time_df_s = Time_df_v1(Start_Time:Stop_Time,:);
        Time_real = Time_df_s - Time_df_s(1);
        

%% Filter
markers_df_s_filt=[];
fc = 1;
fs = 240;
[b,a] = butter(4,fc/(fs/2));
dtm.LeftBeltSpeed_s_filt = filtfilt(b,a,dtm.LeftBeltSpeed_s);
dtm.RightBeltSpeed_s_filt =filtfilt(b,a,dtm.RightBeltSpeed_s);


for i=1:length(markers_df.labels)
    for e=1:3
        markers_df_s_filt.(markers_df.labels{i})(:,e) = filtfilt(b,a,markers_df_s.(markers_df.labels{i})(:,e));
    end
end
       
                        
        %% GET GAIT EVENTS
        HSrefinePre=10;
        HSrefinePost= 5;
        TOminpeakheight=4 ;
        TOminpeakdistance=40;
        BW=0;
        Time = Time_df_s;
        markers4GE = {'RHEE' 'LHEE' 'RANK' 'LANK' 'RTOE' 'LTOE'};
        markers_df.labels(1,:) = {'LASI', 'RASI', 'LPSI', 'RPSI', 'LKNE', 'LTHI', 'LANK', 'LTIB', 'LTOE', 'LHEE', 'RKNE', 'RTHI', 'RANK', 'RTIB', 'RTOE', 'RHEE'};
        llmarkers= markers_df_s_filt;
        
        
%         [RHS,LTO,LHS,RTO,GE] = GaitEvents_allslopes_v7(Time_real,llmarkers,markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW); 
[RHS,LTO,LHS,RTO,GE] = GaitEvents_mocap_v7(Time, llmarkers, markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW,subjs(s), c); 
        
        RHS = RHS';
        LTO = LTO';
        LHS = LHS';
        RTO = RTO';
        
% if strcmp(subjs(s),'SPP16') && strcmp(conds_f(c),'2')
%    RHS(1) = [];
%     GE(:,1)=RHS;
%     GE(:,2)=LTO;
%     GE(:,3)=LHS;
%     GE(:,4)=RTO;
% end
% if strcmp(subjs(s),'SPP16') && strcmp(conds_f(c),'4')
%    RHS(1) = [];
%    RHS = [RHS(1:79); 16250; RHS(80:394)];
%    RTO = [RTO(1:78); 16205; RTO(79:394)];
%     GE(:,1)=RHS;
%     GE(:,2)=LTO;
%     GE(:,3)=LHS;
%     GE(:,4)=RTO;
% end
% if strcmp(subjs(s),'SPP3') && strcmp(conds_f(c),'4')
%    RHS(1) = [];
%     GE(:,1)=RHS;
%     GE(:,2)=LTO;
%     GE(:,3)=LHS;
%     GE(:,4)=RTO; 
% end

if strcmp(subjs(s),'SPP2') && strcmp(conds_f(c),'3')
   RHS(1) = [];
    GE(:,1)=RHS;
    GE(:,2)=LTO;
    GE(:,3)=LHS;
    GE(:,4)=RTO; 
end

        GEgood=GE;
       
%         acceptable = GEgood(:,1)>7200;
%         GEgood = GEgood(acceptable(:),:);
        PPP=[];
        leftheel_GE = [];
        
        PPP=repmat(GEgood(2:end,1),1);
        leftheel_GE = repmat(GEgood(2:end,3),1);
        GEgood(end,:)=[];
        GEgood(:,5)=PPP;
        
        
        %% Convert to conventional coords
        for m = 1:length(markers_df_s.labels)
            eval(['markers_df_c.' markers_df_s.labels{m} ' = convert_coords2conventional(markers_df_s_filt. ' markers_df_s.labels{m} ');']);
        end
        
        for m = 1:length(forces_df.labels)
            eval(['forces_df_c.' forces_df.forces_df_s.labels{m} ' = convert_coords2conventional(forces_df.forces_df_s. ' forces_df.labels{m} ');']);
        end     
        
        %% Gait event check with force plates 
%         figure
%         plot(forces_df_c.FP2For(:,3),'r');
%         hold on
%         plot(forces_df_c.FP1For(:,3),'b');
%         %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k', 'linewidth', 2)
%         %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k')
%         plot(GEgood(:,1), zeros(size(GEgood(:,1))), 'rx', GEgood(:,2), zeros(size(GEgood(:,1))), 'bo', GEgood(:,3), zeros(size(GEgood(:,1))), 'bx', GEgood(:,4), zeros(size(GEgood(:,1))), 'ro', GEgood(:,5), zeros(size(GEgood(:,1))), 'rx')
%         ylimits = ylim(gca);
%         title([subjs{s} ', ' conds{c}])
        %             axis([0 RHS(10) ylimits(1) ylimits(2)]);
% % %         savefig(h,'Incline_gait_ACC.fig')
% % %         close(gcf)
%% analysis
COM_x=(markers_df_c.RASI(:,1)+markers_df_c.LASI(:,1)+markers_df_c.RPSI(:,1)+markers_df_c.LPSI(:,1))/4;
COM.(conds{c}).(subjs{s})(:,1) = COM_x;
COM.(conds{c}).avg(s,1) = mean(COM_x);
COM.(conds{c}).standdev(s,1) = std(COM_x);
COM.(conds{c}).vardev(s,1) = var(COM_x);

COM_all=(markers_df_c.RASI(:,2)+markers_df_c.LASI(:,2)+markers_df_c.RPSI(:,2)+markers_df_c.LPSI(:,2))/4;
COM.(conds{c}).(subjs{s})(:,2) = COM_all;
COM.(conds{c}).avg(s,2) = mean(COM_all);
COM.(conds{c}).standdev(s,2) = std(COM_all);
COM.(conds{c}).vardev(s,2) = var(COM_all);



COM_all_vel=[];
for i = 1:length(Time_real)-1
    COM_all_vel(i) = (COM_all(i+1)-COM_all(i))/(Time_real(i+1)-Time_real(i)) ;
end

COM_all_vel_all{c,s} = COM_all_vel;

COM_all_plus_speed = [];
COM_all_plus_speed_Right = dtm.RightBeltSpeed_s_filt(1:end-1,1) + COM_all_vel';
COM_all_plus_speed_Left = dtm.LeftBeltSpeed_s_filt(1:end-1,1) + COM_all_vel';


      for i = 1:length(GEgood(:,1))
            steplength_speed(i,1) = nanmean(COM_all_plus_speed_Right(GEgood(i,1):GEgood(i,3)));
            steplength_speed(i,2) = nanmean(COM_all_plus_speed_Left(GEgood(i,3):GEgood(i,5)));
      end  
      
      for i = 1:length(GEgood(:,1))
            steplength_time(i,1) = range(Time_df_s(GEgood(i,1):GEgood(i,3)));
            steplength_time(i,2) = range(Time_df_s(GEgood(i,3):GEgood(i,5)));
      end
      
        sl_temp(:,1) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3), GEgood(:,1),steplength_speed(:,1), steplength_time(:,1),0);
        sl_temp(:,2) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,5), GEgood(:,3),steplength_speed(:,2), steplength_time(:,2),0);
        sw_temp(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,3));
        sw_temp(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,5));

        
        steplength_speed_all = [];
        for i = 1:size(steplength_speed,1)
            steplength_speed_all = [steplength_speed_all steplength_speed(i,:)];
        end

        sl_all = [];
        for i = 1:size(sl_temp,1)
            sl_all = [sl_all sl_temp(i,:)];
        end
        steplength_time_all = [];
        for i = 1:size(steplength_time,1)
          steplength_time_all = [steplength_time_all steplength_time(i,:)];
        end
        sw_all = [];
        for i = 1:size(sw_temp,1)
            sw_all = [sw_all sw_temp(i,:)];
        end
        
        ws.ss=size(steplength_speed_all,2);
        ws.(conds{c}).(subjs{s}) = steplength_speed_all; 
%         ws.(conds{c}).([subjs{s} 'stats'])(:,1) = steplength_speed_all(10:440);
%         ws.(conds{c}).([subjs{s} 'statssubj'])(:,1) = repmat(subjs(s),length(steplength_speed_all(10:440)),1);
        ws.(conds{c}).avg(s) = mean(steplength_speed_all);
        ws.(conds{c}).standdev(s) = std(steplength_speed_all);
        ws.(conds{c}).vardev(s) = var(steplength_speed_all);
        ws.(subjs{s}).avg(c) = mean(steplength_speed_all);
        ws.(subjs{s}).standdev(c) = std(steplength_speed_all);
        ws.(subjs{s}).vardev(c) = var(steplength_speed_all);
        ws.(subjs{s}).avg_halfs(c,:) = [mean(steplength_speed_all(1:(ws.ss/2))) mean(steplength_speed_all((ws.ss/2):end))];
        ws.(subjs{s}).standdev_halfs(c,:) = [std(steplength_speed_all(1:(ws.ss/2))) std(steplength_speed_all((ws.ss/2):end))];


        
        sl.ss=size(sl_all,2);        
        sl.(conds{c}).(subjs{s}) = sl_all;
        sl.(conds{c}).avg(s) = mean(sl_all);
        sl.(conds{c}).standdev(s) = std(sl_all);
        sl.(conds{c}).vardev(s) = var(sl_all);
        sl.(subjs{s}).avg(c) = mean(sl_all);
        sl.(subjs{s}).standdev(c) = std(sl_all);
        sl.(subjs{s}).vardev(c) = var(sl_all);
        sl.(subjs{s}).avg_halfs(c,:) = [mean(sl_all(1:(sl.ss/2))) mean(sl_all((sl.ss/2):end))];
        sl.(subjs{s}).standdev_halfs(c,:) = [std(sl_all(1:(sl.ss/2))) std(sl_all((sl.ss/2):end))];


        sw.ss=size(sw_all,2);
        sw.(conds{c}).(subjs{s}) = sw_all; 
        sw.(conds{c}).avg(s) = mean(sw_all);
        sw.(conds{c}).standdev(s) = std(sw_all);
        sw.(conds{c}).vardev(s) = var(sw_all);
        sw.(subjs{s}).avg(c) = mean(sw_all);
        sw.(subjs{s}).standdev(c) = std(sw_all);
        sw.(subjs{s}).vardev(c) = var(sw_all);
        sw.(subjs{s}).avg_halfs(c,:) = [mean(1./sw_all(1:(sw.ss/2))) mean(sw_all((sw.ss/2):end))];
        sw.(subjs{s}).standdev_halfs(c,:) = [std(1./sw_all(1:(sw.ss/2))) std(sw_all((sw.ss/2):end))];
     
        
        sf.ss=size(steplength_time_all,2);
        sf.(conds{c}).(subjs{s}) = 1./steplength_time_all;
        sf.(conds{c}).avg(s) = mean(1./steplength_time_all);
        sf.(conds{c}).standdev(s) = std(1./steplength_time_all);
        sf.(conds{c}).vardev(s) = var(1./steplength_time_all);
        sf.(subjs{s}).avg(c) = mean(1./steplength_time_all);
        sf.(subjs{s}).standdev(c) = std(1./steplength_time_all);
        sf.(subjs{s}).vardev(c) = var(1./steplength_time_all);
        sf.(subjs{s}).avg_halfs(c,:) = [mean(1./steplength_time_all(1:(sf.ss/2))) mean(1./steplength_time_all((sf.ss/2):end))];
        sf.(subjs{s}).standdev_halfs(c,:) = [std(1./steplength_time_all(1:(sf.ss/2))) std(1./steplength_time_all((sf.ss/2):end))];


%% plots 
    
%         figure(3)
%         subplot(2,3,c)
% %         ylim([0.4 0.9])
%         plot(sl.(conds{c}).(subjs{s})); hold on,
% %         ylim([0.4 0.9])
%         title((conds{c}))
%         
%         figure(4)
%         subplot(2,3,c)
%         plot(sf.(conds{c}).(subjs{s})); hold on,
% %         ylim([1.6 2.4])
%         title((conds{c}))
%         
%         figure(5)
%         subplot(2,3,c)
%         plot(sw.(conds{c}).(subjs{s})); hold on,
% %         ylim([0 0.3])
%         title((conds{c}))
%         
%         figure(6)
%         subplot(2,3,c)
%         plot(ws.(conds{c}).(subjs{s})); hold on,
% %         ylim([0 0.3])
%         title((conds{c}))
       
%         figure(7)
%         subplot(2,3,c)
%         plot(dtm.sway_s,'k-'), hold on
%         plot(GEgood(1:4:size(GEgood(:,1)),1),zeros(size(GEgood(1:4:size(GEgood(:,1)),1))),'r.')
% 
%         
        
%         figure(8)
%         subplot(2,3,c)
%         for i=1:2:size(GEgood(:,1))-2
%         plot(1:((GEgood(i+1,1)-GEgood(i,1)+1)),dtm.sway_s(GEgood(i,1):GEgood(i+1,1))), hold on
%         end
%         if c==4 
%             figure(3)
%             subplot(2,3,c)
%             ylabel('step length')
%             figure(4)
%             subplot(2,3,c)
%             ylabel('step frequency')
%             figure(5)
%             subplot(2,3,c)
%             ylabel('step width')
%             figure(6)
%             subplot(2,3,c)
%             ylabel('walkiing speed')
%         end
    end
end
%% Clean workspace
clear a b bins BW COM_all_plus_speed counts Dflow_size e fc HSrefinePost HSrefinePre i istart m NOG Start_Time GE llmarkers P Time_df Time_df_s Time_df_v1 Time_real TOminpeakdistance TOminpeakheight sl_temp Total tf sw_temp Stop_Time stopidx startidx start_values start_value_pref dtm dflow_file dflow_treadmill_file fh;
clear fitplot_actual_steplength fitplot_fitted_steplength fitplot_speed fitplot_stepfitted_minus_actualstep forces_df start_vlue_pref start_values startidx steplength_speed steplength_speed_all steplength_time steplength_time_all Stop_Time stopidx Frame_df Frame_df_s Frame_df_v1 LHS LTO RHS RTO ppp ;      
%% string variables for plots
subjs_s = ["SPP2" "SPP3" "SPP4" "SPP5" "SPP6" "SPP8" "SPP9" "SPP10" "SPP11" "SPP12"];
conds_s= ["no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"];

%% grouping kinematics 
kinem.sl = sl;
kinem.sf = sf;
kinem.sw = sw;
kinem.ws = ws;

save pre_kinem_o pre_kinem_o
%% detrended plot 
m = 1;
variation_steps.stack_p=[];
for i=conds_s    
variation_steps.(i).det_avg = mean(variation_steps.(i).slminusfit);
variation_steps.(i).speedt_avg = mean(variation_steps.(i).speedtrend);
variation_steps.stack_p(m,:) = [mean(variation_steps.(i).slminusfit), mean(variation_steps.(i).speedtrend)];
m = m + 1;
end
figure
bar(variation_steps.stack_p,'stacked')
legend('detrended','speedtrend')
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no pert' 'same mf' 'diff f' 'diff m' 'diff fm'});
ylabel('variance (m^2)')
%% plots
% for m=1:length(subjs)
    for i=1:length(conds)
        figure(100)
        subplot(4,2,1)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i,sl.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sl.(conds{i}).avg),mean(sl.(conds{i}).standdev),'k.'), hold on
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
        plot(i,mean(sl.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i,sl.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sl.(conds{i}).standdev),mean(std(sl.(conds{i}).standdev)),'k.'), hold on
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
        plot(i,mean(sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i,sf.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sf.(conds{i}).avg),mean(sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency avg')
% errorbar([1; 2; 3; 4; 5],mean(sf.(conds{i}).avg),mean(sf.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,4)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sf.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i,sf.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sf.(conds{i}).standdev),mean(std(sf.(conds{i}).standdev)),'k.'), hold on
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
        plot(i,mean(sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i,sw.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sw.(conds{i}).avg),mean(sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')
% errorbar([1; 2; 3; 4; 5],mean(sw.(conds{i}).avg),mean(sw.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,6)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sw.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i,sw.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i,mean(sw.(conds{i}).standdev),mean(std(sw.(conds{i}).standdev)),'k.'), hold on
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
        plot(i,mean(ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
        plot(i,ws.(conds{i}).avg,'b.','MarkerSize',4), hold on
        errorbar(i,mean(ws.(conds{i}).avg),mean(ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')
% errorbar([1; 2; 3; 4; 5],mean(ws.(conds{i}).avg),mean(ws.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(100)
        subplot(4,2,8)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(ws.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        plot(i,ws.(conds{i}).standdev,'b.','MarkerSize',4), hold on
        errorbar(i,mean(ws.(conds{i}).standdev),mean(std(ws.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(100)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')


% for m=1:length(subjs)
%     figure(101)
%     subplot(4,2,2)
% %     title(subjs{m})
% 
%     for i=1:length(conds)
% %         figure(200)
% %         xlim([0 6])
% %         ylim([0.6 0.85])
%         plot([1 2],sl.(subjs{m}).standdev_halfs,'o-'), hold on
%         
%         xlim([0 3])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
% %         plot(i,sl.(conds{i}).avg(m),'o'), hold on
% %         plot(i,sl.(conds{i}).avg(2),'ro'), hold on
% %         legend(conds)
%     end
%     title(subjs{m})
% end