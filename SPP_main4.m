%% NEW CODE 
close all
clearvars
clc

% subjs = {'SPP2' 'SPP3' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11'};
% subjs = {'SPP2' 'SPP3' 'SPP4' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11' 'SPP12'};
% subjs = {'SPP13' 'SPP14' 'SPP15' 'SPP16' 'SPP17' 'SPP19' 'SPP21' 'SPP22'};

subjs = {'SPP23'};

conds_f = {'0' '1' '2' '3' '4'};
% conds_f = {'4'};

conds = {'no_pert' 'same_mf' 'diff_f' 'diff_m' 'diff_fm'};
% conds = {'no_pert'};


fs = 240; %hz, resamp tm to df rate
proj = 'Z:\SPP\subjects\';
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
    istart = find(dtm.sway_v1 < -0.005,1,'first');
    Stop_Time = find(dtm.stride_v1 < 480,1,'last');
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
        
if strcmp(subjs(s),'SPP16') && strcmp(conds_f(c),'2')
   RHS(1) = [];
    GE(:,1)=RHS;
    GE(:,2)=LTO;
    GE(:,3)=LHS;
    GE(:,4)=RTO;
end
if strcmp(subjs(s),'SPP16') && strcmp(conds_f(c),'4')
   RHS(1) = [];
   RHS = [RHS(1:79); 16250; RHS(80:394)];
   RTO = [RTO(1:78); 16205; RTO(79:394)];
    GE(:,1)=RHS;
    GE(:,2)=LTO;
    GE(:,3)=LHS;
    GE(:,4)=RTO;
end
if strcmp(subjs(s),'SPP3') && strcmp(conds_f(c),'4')
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

% detrended analysis      

if strcmp((subjs{s}),('SPP13')) && strcmp((conds{s}),('no_pert'))
        NOG = (length(steplength_speed_all)) - 60;
        
        
        [P,fh] = fitsinglemodelprocess_sl(sl.(conds{c}).(subjs{s})(NOG:end),ws.(conds{c}).(subjs{s})(NOG:end));
        p.(conds{c})(s,:) = P;
        
        fitplot_fitted_steplength.(conds{c}).(subjs{s}) = fh(ws.(conds{c}).(subjs{s})(NOG:end),P);
        fitplot_stepfitted_minus_actualstep.(conds{c}).(subjs{s}) = (sl.(conds{c}).(subjs{s})(NOG:end) - fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        fitplot_speed.(conds{c}).(subjs{s}) = ws.(conds{c}).(subjs{s})(NOG:end);
        fitplot_actual_steplength.(conds{c}).(subjs{s}) = sl.(conds{c}).(subjs{s})(NOG:end);
        
        variation_steps.(conds{c}).slminusfit(s,:) = var(sl.(conds{c}).(subjs{s})(NOG:end) - fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        variation_steps.(conds{c}).speedtrend(s,:) = var(fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        variation_steps.(conds{c}).totalvar(s,:) = var(sl.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).totalstd(s,:) = std(sl.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).speedvar(s,:) = var(ws.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).speedmean(s,:) = mean(ws.(conds{c}).(subjs{s})(NOG:end));
else
        NOG = (length(steplength_speed_all)) - 400;
        
        
        [P,fh] = fitsinglemodelprocess_sl(sl.(conds{c}).(subjs{s})(NOG:end),ws.(conds{c}).(subjs{s})(NOG:end));
        p.(conds{c})(s,:) = P;
        
        fitplot_fitted_steplength.(conds{c}).(subjs{s}) = fh(ws.(conds{c}).(subjs{s})(NOG:end),P);
        fitplot_stepfitted_minus_actualstep.(conds{c}).(subjs{s}) = (sl.(conds{c}).(subjs{s})(NOG:end) - fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        fitplot_speed.(conds{c}).(subjs{s}) = ws.(conds{c}).(subjs{s})(NOG:end);
        fitplot_actual_steplength.(conds{c}).(subjs{s}) = sl.(conds{c}).(subjs{s})(NOG:end);
        
        variation_steps.(conds{c}).slminusfit(s,:) = var(sl.(conds{c}).(subjs{s})(NOG:end) - fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        variation_steps.(conds{c}).speedtrend(s,:) = var(fh(ws.(conds{c}).(subjs{s})(NOG:end),P));
        variation_steps.(conds{c}).totalvar(s,:) = var(sl.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).totalstd(s,:) = std(sl.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).speedvar(s,:) = var(ws.(conds{c}).(subjs{s})(NOG:end));
        variation_steps.(conds{c}).speedmean(s,:) = mean(ws.(conds{c}).(subjs{s})(NOG:end));
end

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
save main_kinem_o kinem_o
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
%% COM plot 
n=1;
for i=conds_s
    for m=subjs_s
    figure(200)
    subplot(2,3,n)
    title(conds_s(n))
    plot(COM.(i).(m)); hold on
    end
n=n+1;
end   
%% more plots 
figure(200)
% set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
% title('step length std halfs')
legend(conds)

for m=1:length(subjs)
    figure(201)
%     subplot(2,3,m)
%     title(subjs{m})

%     for i=1:length(conds)
        plot(sl.(subjs{m}).standdev_halfs(:,2),'o-'), hold on
        xlim([0 6])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
%     end
%     title(subjs{m})
end
figure(201)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step length std 2nd half')
legend(subjs)

for m=1:length(subjs)
    figure(202)
%     subplot(2,3,m)
%     title(subjs{m})

%     for i=1:length(conds)
        plot(sl.(subjs{m}).avg_halfs(:,2),'o-'), hold on
        xlim([0 6])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
%     end
%     title(subjs{m})
end
figure(202)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step length avg 2nd half')
legend(subjs)

for m=1:length(subjs)
    figure(203)
%     subplot(2,3,m)
%     title(subjs{m})

%     for i=1:length(conds)
        plot(sw.(subjs{m}).standdev_halfs(:,2),'o-'), hold on
        xlim([0 6])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
%     end
%     title(subjs{m})
end
figure(203)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step width std 2nd half')
legend(subjs)

for m=1:length(subjs)
    figure(203)
%     subplot(2,3,m)
%     title(subjs{m})

%     for i=1:length(conds)
        plot(sw.(subjs{m}).standdev_halfs(:,2),'o-'), hold on
        xlim([0 6])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
%     end
%     title(subjs{m})
end
figure(203)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step width std 2nd half')
legend(subjs)

for m=1:length(subjs)
    figure(203)
%     subplot(2,3,m)
%     title(subjs{m})

%     for i=1:length(conds)
        plot(sw.(subjs{m}).standdev_halfs(:,2),'o-'), hold on
        xlim([0 6])
%         set(gca, 'XTick', [1,2],'XTickLabel',{'first' 'second'});
%     end
%     title(subjs{m})
end
figure(203)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step width std 2nd half')
legend(subjs)


% for m=1:length(subjs)
%     for i=1:length(conds)
%         figure(101)
%         xlim([0 6])
% %         ylim([0.6 0.85])
%         plot(i,sl.(conds{i}).standdev(m),'o'), hold on
% %         plot(i,sl.(conds{i}).avg(2),'ro'), hold on
%         legend(subjs)
%     end
% end
% figure(101)
% set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
% title('step length std')

% for m=1:length(subjs)
% %     for i=1:length(conds)
%         figure(102)
%         xlim([0 6])
%         ylim([1.8 2.2])
%         plot(sf.(subjs{m}).avg,'o-'), hold on
% 
% %         plot(i,sf.(conds{i}).avg(m),'bo'), hold on
% %         plot(i,sf.(conds{i}).avg(2),'ro'), hold on
%         legend(subjs)
% %     end
% end
% figure(102)
% set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
% title('step frequency')

for m=1:length(subjs)
%     for i=1:length(conds)
        figure(103)
        xlim([0 6])
%         ylim([1.8 2.2])
         plot(sl.(subjs{m}).avg,'o-'), hold on

%         plot(i,sf.(conds{i}).standdev(m),'bo'), hold on
%         plot(i,sf.(conds{i}).avg(2),'ro'), hold on
        legend(subjs)
%     end
end
figure(103)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step length avg')


% for m=1:length(subjs)
%     for i=1:length(conds)
%         figure(104)
%         xlim([0 6])
%         ylim([0.04 0.2])
%         plot(i,sw.(conds{i}).avg(m),'o'), hold on
% %         plot(i,sw.(conds{i}).avg(2),'ro'), hold on
%         legend(subjs)
%     end
% end
% figure(104)
% set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
% title('step width')

for m=1:length(subjs)
%     for i=1:length(conds)
        figure(105)
        xlim([0 6])
%         ylim([0.04 0.2])
%         plot(i,sw.(conds{i}).avg(m),'o'), hold on
         plot(sw.(subjs{m}).avg,'o-'), hold on
%          plot(i, sw.(conds{i}).avg(2),'o-','Color',[0.8500 0.3250 0.0980]), hold on
%          plot(i, sw.(conds{i}).avg(3),'o-','Color',[0.9290 0.6940 0.1250]), hold on
%          plot(i, sw.(conds{i}).avg(4),'o-','Color',[0.4940 0.1840 0.5560]), hold on
%         plot(i,sw.(conds{i}).avg(2),'ro'), hold on
        legend(subjs)
%     end
end
figure(105)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step width avg')
 
%  for m=1:length(subjs)
     for i=1:length(conds)
         figure(106)
         xlim([0 6])
%          ylim([1.00 1.5])
         plot(i, ws.(conds{i}).avg(1),'bo'), hold on
         plot(i, ws.(conds{i}).avg(2),'ro'), hold on
         plot(i, ws.(conds{i}).avg(3),'ko'), hold on
         plot(i, ws.(conds{i}).avg(4),'go'), hold on
         legend(subjs)
     end
%  end     
figure(106)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('Walking speed')

for m=1:length(subjs)
    
         figure(108)
         xlim([0 6])
%          ylim([1.00 1.5])
         plot(sw.(subjs{m}).standdev,'-o'), hold on
%          plot(i, ws.(conds{i}).standdev(2),'ro'), hold on
%          plot(i, ws.(conds{i}).standdev(3),'ko'), hold on
%          plot(i, ws.(conds{i}).standdev(4),'go'), hold on
         legend(subjs)
 
 end     
figure(108)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('step width std')

for m=1:length(subjs)
    
         figure(207)
         xlim([0 6])
         ylim([1.00 1.5])
         plot(ws.(subjs{m}).avg,'-o'), hold on
%          plot(i, ws.(conds{i}).standdev(2),'ro'), hold on
%          plot(i, ws.(conds{i}).standdev(3),'ko'), hold on
%          plot(i, ws.(conds{i}).standdev(4),'go'), hold on
         legend(subjs)
 
 end     
figure(207)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('Walking speed avg')

for m=1:length(subjs)
    
         figure(109)
         xlim([0 6])
%          ylim([1.00 1.5])
         plot(sf.(subjs{m}).standdev,'-o'), hold on
%          plot(i, ws.(conds{i}).standdev(2),'ro'), hold on
%          plot(i, ws.(conds{i}).standdev(3),'ko'), hold on
%          plot(i, ws.(conds{i}).standdev(4),'go'), hold on
         legend(subjs)
 
 end     
figure(109)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('Step freq std')

for m=1:length(subjs)
    
         figure(110)
         xlim([0 6])
%          ylim([1.00 1.5])
         plot(sl.(subjs{m}).standdev,'-o'), hold on
%          plot(i, ws.(conds{i}).standdev(2),'ro'), hold on
%          plot(i, ws.(conds{i}).standdev(3),'ko'), hold on
%          plot(i, ws.(conds{i}).standdev(4),'go'), hold on
         legend(subjs)
 
 end     
figure(110)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'samefm' 'f' 'm' 'fm'});
title('Step length std')

% tf_filt.speed_afterp.avg(1) = [];
% figure
% plot([1,2,3,4],tf_filt.speed_afterp.avg,'.','MarkerSize',30)
% xlim([0 5])
% set(gca, 'XTick', [1,2,3,4],'XTickLabel',{'nofm' 'f' 'm' 'fm'});
% 
% % 
% % 
% figure(1)
% legend(conds)
% 
% figure(2)
% legend(conds{2:5})
% 
% % 
% % figure(1)
% % plot(tf_filt.(conds{c})); hold on,
% % plot(tf.(conds{c}).SpeedActual2);
% 
% 
%% github changes 
% Stats 
% ws.(conds{c}).(subjs{s})
% repmat(subjs(s),length(steplength_speed_all),1)
for m=1:length(subjs) 
ws.([subjs{m} 'alldata']) = [mean(ws.(conds{1}).([subjs{m} 'stats'])(:,1)); ws.(conds{2}).([subjs{m} 'stats'])(:,1); ws.(conds{3}).([subjs{m} 'stats'])(:,1); ws.(conds{4}).([subjs{m} 'stats'])(:,1); ws.(conds{5}).([subjs{m} 'stats'])(:,1)];   
end
for m=1 
ws.([subjs{m} 'alldatasubjs']) = [ws.(conds{1}).([subjs{m} 'stats'])(:,1); ws.(conds{2}).([subjs{m} 'stats'])(:,1); ws.(conds{3}).([subjs{m} 'stats'])(:,1); ws.(conds{4}).([subjs{m} 'stats'])(:,1); ws.(conds{5}).([subjs{m} 'stats'])(:,1)];   
end

for m=1:length(subjs) 
ws.(['alldataspss'])(m,:) = [mean(ws.(conds{1}).([subjs{m} 'stats'])(:,1)), mean(ws.(conds{2}).([subjs{m} 'stats'])(:,1)), mean(ws.(conds{3}).([subjs{m} 'stats'])(:,1)), mean(ws.(conds{4}).([subjs{m} 'stats'])(:,1)), mean(ws.(conds{5}).([subjs{m} 'stats'])(:,1))];   
end
for m=1 
ws.([subjs{m} 'alldatasubjs']) = [ws.(conds{1}).([subjs{m} 'stats'])(:,1); ws.(conds{2}).([subjs{m} 'stats'])(:,1); ws.(conds{3}).([subjs{m} 'stats'])(:,1); ws.(conds{4}).([subjs{m} 'stats'])(:,1); ws.(conds{5}).([subjs{m} 'stats'])(:,1)];   
end
ws.('subjsalldataconds') = [repmat(conds(1),length(steplength_speed_all(10:475)),1);repmat(conds(2),length(steplength_speed_all(10:475)),1);repmat(conds(3),length(steplength_speed_all(10:475)),1);repmat(conds(4),length(steplength_speed_all(10:475)),1);repmat(conds(5),length(steplength_speed_all(10:475)),1)];

ws.('spss') = [repmat(1,length(steplength_speed_all(10:475)),1);repmat(2,length(steplength_speed_all(10:475)),1);repmat(3,length(steplength_speed_all(10:475)),1);repmat(4,length(steplength_speed_all(10:475)),1);repmat(5,length(steplength_speed_all(10:475)),1);repmat(6,length(steplength_speed_all(10:475)),1);repmat(7,length(steplength_speed_all(10:475)),1);repmat(8,length(steplength_speed_all(10:475)),1)];


for i=1:length(conds) 
ws.(conds{i}).alldata = [ws.(conds{i}).([subjs{1} 'stats'])(:,1); ws.(conds{i}).([subjs{2} 'stats'])(:,1); ws.(conds{i}).([subjs{3} 'stats'])(:,1); ws.(conds{i}).([subjs{4} 'stats'])(:,1); ws.(conds{i}).([subjs{5} 'stats'])(:,1); ws.(conds{i}).([subjs{6} 'stats'])(:,1); ws.(conds{i}).([subjs{7} 'stats'])(:,1); ws.(conds{i}).([subjs{8} 'stats'])(:,1)];   
end

for i=1:length(conds) 
ws.(conds{i}).allsubjsdata = [ws.(conds{i}).([subjs{1} 'statssubj'])(:,1); ws.(conds{i}).([subjs{2} 'statssubj'])(:,1); ws.(conds{i}).([subjs{3} 'statssubj'])(:,1); ws.(conds{i}).([subjs{4} 'statssubj'])(:,1); ws.(conds{i}).([subjs{5} 'statssubj'])(:,1); ws.(conds{i}).([subjs{6} 'statssubj'])(:,1); ws.(conds{i}).([subjs{7} 'statssubj'])(:,1); ws.(conds{i}).([subjs{8} 'statssubj'])(:,1)];   
end

stats_names = {'speed' 'speedvar' 'stridelength_detrended' 'stridelength_totalvar' 'stridelength_speedtrend' 'stridelength_slopes_detrended' 'stridewidth_slopes_detrended' 'stridewidth_speedtrend' 'stridewidth_detrended' 'stridewidth_totalvar'};
for i=1:length(stats_names)
    for e=1:5
        stats.(stats_names{i}).(predictability{e})=table('Size',[s 3],'VariableTypes',repmat("double",1,3),'VariableNames',["Condition_050" "Condition_100" "Condition_150"]);
    end
end

new=[1,11,21];
for m=1:3
    for i=1:s
        for e=1:3       
            stats.speed.(slopes{m}){i,e}=stride_variation_speedmean_scatter_plot.(slopes{m})(i,e); 
            stats.speedvar.(slopes{m}){i,e}=stride_variation_speedvar_scatter_plot.(slopes{m})(i,e);
            stats.stridelength_detrended.(slopes{m}){i,e}=stridelength_variation_detrended_scatter_plot.(slopes{m})(i,e); 
            stats.stridelength_totalvar.(slopes{m}){i,e}=stridelength_variation_total_scatter_plot.(slopes{m})(i,e);
            stats.stridelength_speedtrend.(slopes{m}){i,e}=stridelength_variation_speedtrend_scatter_plot.(slopes{m})(i,e);
            stats.stridewidth_speedtrend.(slopes{m}){i,e}=stridewidth_variation_speedtrend_scatter_plot.(slopes{m})(i,e);
            stats.stridewidth_totalvar.(slopes{m}){i,e}=stridewidth_variation_total_scatter_plot.(slopes{m})(i,e);
            stats.stridewidth_detrended.(slopes{m}){i,e}=stridewidth_variation_detrended_scatter_plot.(slopes{m})(i,e);
        end
    end
end 


for i=1:length(stats_names)
    for e=1:3
        stats.(stats_names{i}).rm.(slopes{e})=fitrm([stats.(stats_names{i}).(slopes{e})],'Condition_050-Condition_150~1');
        stats.(stats_names{i}).rmResult.(slopes{e})=ranova(stats.(stats_names{i}).rm.(slopes{e}));
        stats.(stats_names{i}).pairtest.(slopes{e})=multCOM_allpare(stats.(stats_names{i}).rm.(slopes{e}),'Time','COM_allparisonType','tukey-kramer');
    end
end

% COM_allpare slopes 

stats.speed_allslopes.slopes_means=table(ws.(conds{1}).allsubjsdata,ws.(conds{1}).alldata,ws.(conds{2}).alldata,ws.(conds{3}).alldata,ws.(conds{4}).alldata,ws.(conds{5}).alldata,'VariableNames',["subjs" "no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"]);
% stats.speed_allslopes.slopes_means=table(ws.('subjsalldataconds'),ws.([subjs{1} 'alldata']),ws.([subjs{2} 'alldata']),ws.([subjs{3} 'alldata']),ws.([subjs{4} 'alldata']),ws.([subjs{5} 'alldata']),ws.([subjs{6} 'alldata']),ws.([subjs{7} 'alldata']),ws.([subjs{8} 'alldata']),'VariableNames',{'conditions','subject1','subject2','subject3','subject4','subject5','subject6','subject7','subject8'});
t=stats.speed_allslopes.slopes_means;
YourArray = table2array(stats.speed_allslopes.slopes_means);
YourNewTable = array2table(YourArray.');
YourNewTable.Properties.VariableNames = conditions';

Xc = table2cell(stats.speed_allslopes.slopes_means);
Xt = cell2table(Xc','RowNames',X.Properties.VariableNames,'VariableNames',X.Properties.RowNames);
% stats.speed_allslopes.slopes_means = table('Size',[s 6],'VariableTypes',repmat("double",1,6),'VariableNames',["subjs" "no_pert" "same_mf" "diff_f" "diff_m" "diff_fm"]);
% stats.speed_allslopes.slopes_means{:,1} = ['SPP2' 'SPP3' 'SPP5' 'SPP6' 'SPP8' 'SPP9' 'SPP10' 'SPP11']';
Meas = table([1 2 3 4 5 6 7 8]','VariableNames',{'subjss'});
rm = fitrm(stats.speed_allslopes.slopes_means,'subject1-subject8~conditions','WithinDesign',Meas);
ranovatbl = ranova(rm);
for e=1:5 
stats.speed_allslopes.slopes_means{:,e+1} = ws.(conds{e}).avg';
end 
% Meas = table([1 2 3 4 5 6 7 8]','VariableNames',{'conditions'});
for e=1:5
        stats.speed_allslopes.rm.slopes_means=fitrm(stats.speed_allslopes.slopes_means,'no_pert-diff_fm~1','WithinDesign');
        stats.speed_allslopes.rmResult.slopes_means=ranova(stats.speed_allslopes.rm.slopes_means);
        stats.speed_allslopes.pairtest.slopes_means=multCOM_allpare(rm,'conditions','COM_allparisonType','tukey-kramer');
end

%% NACOB fig 
color = jet(length(subjs)); 


 for i=1:length(conds)
        figure(99)
        subplot(4,2,1)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sl.(conds{i}).avg),'k.','MarkerSize',20), hold on
        for m=1:length(subjs)
            plot(i,sl.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sl.(conds{i}).avg),mean(sl.(conds{i}).standdev),'k.'), hold on
 end 
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length avg')
% errorbar([1; 2; 3; 4; 5],mean(sl.(conds{i}).avg),mean(sl.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(99)
        subplot(4,2,2)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sl.(conds{i}).standdev),'k.','MarkerSize',20), hold on
  for m=1:length(subjs)
            plot(i,sl.(conds{i}).standdev(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
  end
  errorbar(i,mean(sl.(conds{i}).standdev),mean(std(sl.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step length std')

    for i=1:length(conds)
        figure(99)
        subplot(4,2,3)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sf.(conds{i}).avg),'k.','MarkerSize',20), hold on
        for m=1:length(subjs)
            plot(i,sf.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sf.(conds{i}).avg),mean(sf.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency avg')
% errorbar([1; 2; 3; 4; 5],mean(sf.(conds{i}).avg),mean(sf.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(99)
        subplot(4,2,4)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sf.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        for m=1:length(subjs)
            plot(i,sf.(conds{i}).standdev(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sf.(conds{i}).standdev),mean(std(sf.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step frequency std')

    for i=1:length(conds)
        figure(99)
        subplot(4,2,5)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sw.(conds{i}).avg),'k.','MarkerSize',20), hold on
        for m=1:length(subjs)
            plot(i,sw.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sw.(conds{i}).avg),mean(sw.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width avg')
% errorbar([1; 2; 3; 4; 5],mean(sw.(conds{i}).avg),mean(sw.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(99)
        subplot(4,2,6)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(sw.(conds{i}).standdev),'k.','MarkerSize',20), hold on
         for m=1:length(subjs)
            plot(i,sw.(conds{i}).standdev(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(sw.(conds{i}).standdev),mean(std(sw.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('step width std')

for i=1:length(conds)
        figure(99)
        subplot(4,2,7)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(ws.(conds{i}).avg),'k.','MarkerSize',20), hold on
         for m=1:length(subjs)
            plot(i,ws.(conds{i}).avg(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(ws.(conds{i}).avg),mean(ws.(conds{i}).standdev),'k.'), hold on
    end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed avg')
% errorbar([1; 2; 3; 4; 5],mean(ws.(conds{i}).avg),mean(ws.(conds{i}).standdev),'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
% errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
% [i-0.1,i,i+0.1]
for i=1:length(conds)
        figure(99)
        subplot(4,2,8)
        xlim([0 6])
%         ylim([0.6 0.85])
        plot(i,mean(ws.(conds{i}).standdev),'k.','MarkerSize',20), hold on
        for m=1:length(subjs)
            plot(i,ws.(conds{i}).standdev(m),'.','MarkerSize',10,'Color',color(m,:)), hold on
        end
        errorbar(i,mean(ws.(conds{i}).standdev),mean(std(ws.(conds{i}).standdev)),'k.'), hold on
end
% end
figure(99)
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'nop' 'same t&m' 'dif time' 'dif mag' 'dif t&m'});
title('walking speed std')

%% abstract
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

%% detrended plot 

figure(30) 
for m=1
for i=1:3
subplot(4,4,i)
yyaxis left 
plot(fitplot_actual_steplength.(conds{i}).(subjs{m}),'-','Color','k','LineWidth',1)
ylim([0.4 0.8])
ylabel('Actual Step length (m)')
yyaxis right
plot(fitplot_speed.(conds{i}).(subjs{m}),'-','LineWidth',1) 
ylabel('Walking speed (m/s)')
xlabel('Steps')
% ylim([0.5 1.2])
% xlim([0 325])
title((conds{i}))

% figure
subplot(4,4,4+i)
plot(0.7:0.05:1.3,fh(0.7:0.05:1.3,p.(conds{i})(m,:)),'--'),hold on
plot(fitplot_speed.(conds{i}).(subjs{m}),fitplot_actual_steplength.(conds{i}).(subjs{m}),'.','Color','k','MarkerSize',12),hold on
plot(fitplot_speed.(conds{i}).(subjs{m}),fitplot_actual_steplength.(conds{i}).(subjs{m}),'-','Color','k'), hold on
[counts,bins] = hist(fitplot_actual_steplength.(conds{i}).(subjs{m}));
counts = (0.001 .* counts) + 0.7;
barh(bins,counts,'histc')
if i==3
legend('Fit','Actual Step length')
end
ylabel('Length (m)')
xlabel('Speed (m/s)')
% ylim([0.4 0.7])
% xlim([0.7 1.3])

% figure
subplot(4,4,8+i)
plot(fitplot_actual_steplength.(conds{i}).(subjs{m}),'.','Color','k'),hold on
plot(fitplot_fitted_steplength.(conds{i}).(subjs{m}),'O','Color','g','LineWidth',2), hold on  

if i==3
legend('Actual Step length','Speed-trend')
end
% xlim([0 325])
ylabel('Step Length (m)')
xlabel('Steps')
% ylim([0.4 0.7])


% figure
subplot(4,4,12+i)
plot(fitplot_stepfitted_minus_actualstep.(conds{i}).(subjs{m}),'.','Color','g','MarkerSize',12), hold on
% xlim([0 325])
ylabel('Detrended')
xlabel('Steps')
% ylim([-0.075 0.075])

end
end

% for i=1:3
% % figure(30) 
% % figure1 = figure('WindowState','maximized');
% box(subplot(4,4,8+i),'on');
% txt_on_value2=num2str(variation_steps.(conds{i}).speedtrend(m,:));
% merge_text2={'Speed-trend variance',txt_on_value2};
% annotation(figure(30),'textbox',...
%     [0.004125+(i*0.2) 0.330456905503634 0.0495208333333333 0.0218068535825545],...
%     'String',merge_text2,...
%     'FitBoxToText','on');
% end
% 
% for i=1:3
%     
% 
% % figure(30) 
% % figure1 = figure('WindowState','maximized');
% box(subplot(4,4,12+i),'on');
% txt_on_value=num2str(variation_steps.(conds{i}).slminusfit(m,:));
% merge_text={'Detrended variance',txt_on_value};
% annotation(figure(30),'textbox',[0.004125+(i*0.2) 0.12456905503634 0.0495208333333333 0.0218068535825545],...
%     'String',merge_text,...
%     'FitBoxToText','on');
% end

%% VARIANCE PLOTS
% detrended 
for i=[1 2 3 4 5]
        error_length_detrended(i)=std(variation_steps.(conds{i}).slminusfit(:));
end
figure 
subplot(1,3,2)
stridelength_variation_detrended_bar_plot = [mean(variation_steps.no_pert.slminusfit(:)) mean(variation_steps.same_mf.slminusfit(:)) mean(variation_steps.diff_f.slminusfit(:)) mean(variation_steps.diff_m.slminusfit(:)) mean(variation_steps.diff_fm.slminusfit(:))];
bar(stridelength_variation_detrended_bar_plot);
set(gca, 'XTick', [1,2,3,4,5],'XTickLabel',{'no_p' 'sameTM' 'T_diff' 'M_diff' 'TM_diff'});
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_detrended_bar_plot,error_length_detrended,'.');
for i=1:s
stridelength_variation_detrended_scatter_plot.mean = [mean(variation_steps.no_pert.slminusfit(i,:),2),mean(variation_steps.same_mf.slminusfit(i,:),2),mean(variation_steps.diff_f.slminusfit(i,:),2),mean(variation_steps.diff_m.slminusfit(i,:),2),mean(variation_steps.diff_fm.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_detrended_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_detrended_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_detrended_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
% for e=1:5
% stridelength_variation_detrended_scatter_plot.(slopes{e})(i,:) = stridelength_variation_detrended_scatter_plot.mean(e,:);
% end
end 
ylabel('Step Length Variability -- Detrended')
ylim([0 0.002])
% Total Var 
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_totalvar(1,i-6)=std(variation_steps.(conds{i}).totalvar(:));
    elseif (1<=i&&i<=3)
        error_length_totalvar(2,i)=std(variation_steps.(conds{i}).totalvar(:));
    elseif (4<=i&&i<=6)
        error_length_totalvar(3,i-3)=std(variation_steps.(conds{i}).totalvar(:));
    end 
end
subplot(1,3,3) 
stridelength_variation_total_bar_plot = [mean(variation_steps.decline_050.totalvar(:)),mean(variation_steps.decline_100.totalvar(:)),mean(variation_steps.decline_150.totalvar(:));mean(variation_steps.level_050.totalvar(:)),mean(variation_steps.level_100.totalvar(:)),mean(variation_steps.level_150.totalvar(:));mean(variation_steps.incline_050.totalvar(:)),mean(variation_steps.incline_100.totalvar(:)),mean(variation_steps.incline_150.totalvar(:))];
bar(stridelength_variation_total_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_total_bar_plot,error_length_totalvar,'.');
for i=1:s
stridelength_variation_total_scatter_plot.mean = [mean(variation_steps.decline_050.totalvar(i,:),2),mean(variation_steps.decline_100.totalvar(i,:),2),mean(variation_steps.decline_150.totalvar(i,:),2);mean(variation_steps.level_050.totalvar(i,:),2),mean(variation_steps.level_100.totalvar(i,:),2),mean(variation_steps.level_150.totalvar(i,:),2);mean(variation_steps.incline_050.totalvar(i,:),2),mean(variation_steps.incline_100.totalvar(i,:),2),mean(variation_steps.incline_150.totalvar(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_total_scatter_plot.(slopes{e})(i,:) = stridelength_variation_total_scatter_plot.mean(e,:);
end
end 
ylabel('Step Length Variability -- Total')
ylim([0 0.004])
% Speedtrend
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_speedtrend(1,i-6)=std(variation_steps.(conds{i}).speedtrend(:));
    elseif (1<=i&&i<=3)
        error_length_speedtrend(2,i)=std(variation_steps.(conds{i}).speedtrend(:));
    elseif (4<=i&&i<=6)
        error_length_speedtrend(3,i-3)=std(variation_steps.(conds{i}).speedtrend(:));
    end 
end
subplot(1,3,1) 
stridelength_variation_speedtrend_bar_plot = [mean(variation_steps.decline_050.speedtrend(:)),mean(variation_steps.decline_100.speedtrend(:)),mean(variation_steps.decline_150.speedtrend(:));mean(variation_steps.level_050.speedtrend(:)),mean(variation_steps.level_100.speedtrend(:)),mean(variation_steps.level_150.speedtrend(:));mean(variation_steps.incline_050.speedtrend(:)),mean(variation_steps.incline_100.speedtrend(:)),mean(variation_steps.incline_150.speedtrend(:))];
bar(stridelength_variation_speedtrend_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_speedtrend_bar_plot,error_length_speedtrend,'.');
for i=1:s
stridelength_variation_speedtrend_scatter_plot.mean = [mean(variation_steps.decline_050.speedtrend(i,:),2),mean(variation_steps.decline_100.speedtrend(i,:),2),mean(variation_steps.decline_150.speedtrend(i,:),2);mean(variation_steps.level_050.speedtrend(i,:),2),mean(variation_steps.level_100.speedtrend(i,:),2),mean(variation_steps.level_150.speedtrend(i,:),2);mean(variation_steps.incline_050.speedtrend(i,:),2),mean(variation_steps.incline_100.speedtrend(i,:),2),mean(variation_steps.incline_150.speedtrend(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_speedtrend_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_speedtrend_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_speedtrend_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_speedtrend_scatter_plot.(slopes{e})(i,:) = stridelength_variation_speedtrend_scatter_plot.mean(e,:);
end
end 
ylabel('Step Length Variability -- Speedtrend')
ylim([0 0.002])


%% STACKED VAR
% for i=[7 8 9 1 2 3 4 5 6]
%     if (7<=i&&i<=9)
%         error_length_detrended(1,i-6)=std(variation_steps.(conds{i}).slminusfit(:));
%     elseif (1<=i&&i<=3)
%         error_length_detrended(2,i)=std(variation_steps.(conds{i}).slminusfit(:));
%     elseif (4<=i&&i<=6)
%         error_length_detrended(3,i-3)=std(variation_steps.(conds{i}).slminusfit(:));
%     end 
% end
% stridelength_stacked=[];
% for i=[1 4 7 2 5 8 3 6 9]
% e=size(stridelength_stacked,1)+1;
% stridelength_stacked(e,:)=[stridelength_variation_detrended_bar_plot(i),stridelength_variation_speedtrend_bar_plot(i)];
% stridelength_stacked(e+1,:)=[stridelength_variation_total_bar_plot(i),0]; 
% end 



figure;
h=bar([0.80 0.90 1.15 1.25 1.50 1.60],stridelength_stacked(1:6,:),'stacked'); hold on
colors_stacks = jet(size(h,2));
colors_stacks = repelem(colors_stacks,size(h,1),1); 
colors_stacks = mat2cell(colors_stacks,ones(size(colors_stacks,1),1),3);
set(h,{'FaceColor'},colors_stacks)
h=bar([2.60 2.70 2.95 3.05 3.30 3.40],stridelength_stacked(7:12,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
hold on
h=bar([4.40 4.50 4.75 4.85 5.10 5.20],stridelength_stacked(13:18,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
set(gca, 'XTick', [1.2,3,4.8],'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
for i=1:s
stridelength_variation_total_scatter_plot.mean = [mean(variation_steps.decline_050.totalvar(i,:),2),mean(variation_steps.decline_100.totalvar(i,:),2),mean(variation_steps.decline_150.totalvar(i,:),2);mean(variation_steps.level_050.totalvar(i,:),2),mean(variation_steps.level_100.totalvar(i,:),2),mean(variation_steps.level_150.totalvar(i,:),2);mean(variation_steps.incline_050.totalvar(i,:),2),mean(variation_steps.incline_100.totalvar(i,:),2),mean(variation_steps.incline_150.totalvar(i,:),2)];
% plot([0.90 1.25 1.60],stridelength_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([2.70 3.05 3.40],stridelength_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([4.50 4.85 5.20],stridelength_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_total_scatter_plot.(slopes{e})(i,:) = stridelength_variation_total_scatter_plot.mean(e,:);
end
end 
hold on 
errorbar([0.90 1.25 1.60; 2.70 3.05 3.40; 4.50 4.85 5.20],stridelength_variation_total_bar_plot,error_length_totalvar,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridelength_variation_detrended_bar_plot,error_length_detrended,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridelength_variation_detrended_bar_plot + stridelength_variation_speedtrend_bar_plot),error_length_speedtrend,'.');
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Length Variability')


%% DETRENDED VS SLOPE 
figure 
stridelength_variation_detrended_vs_slopes = [mean(variation_steps.decline_050.slminusfit(:)),mean(variation_steps.level_050.slminusfit(:)),mean(variation_steps.incline_050.slminusfit(:));mean(variation_steps.decline_100.slminusfit(:)),mean(variation_steps.level_100.slminusfit(:)),mean(variation_steps.incline_100.slminusfit(:));mean(variation_steps.decline_150.slminusfit(:)),mean(variation_steps.level_150.slminusfit(:)),mean(variation_steps.incline_150.slminusfit(:))];
% bar(stridelength_variation_detrended_vs_slopes)
% set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});

for i=1:s
stridelength_variation_detrended_vs_slope_scatter_plot.mean = [mean(variation_steps.decline_050.slminusfit(i,:),2),mean(variation_steps.level_050.slminusfit(i,:),2),mean(variation_steps.incline_050.slminusfit(i,:),2);mean(variation_steps.decline_100.slminusfit(i,:),2),mean(variation_steps.level_100.slminusfit(i,:),2),mean(variation_steps.incline_100.slminusfit(i,:),2);mean(variation_steps.decline_150.slminusfit(i,:),2),mean(variation_steps.level_150.slminusfit(i,:),2),mean(variation_steps.incline_150.slminusfit(i,:),2)];
plot([0.78 1 1.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',[0,0,0]+0.8); hold on 
plot([1.78 2 2.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',[0,0,0]+0.8)
plot([2.78 3 3.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',[0,0,0]+0.8)
for e=1:3
stridelength_variation_detrended_vs_slope_scatter_plot.(slopes{e})(i,:) = stridelength_variation_detrended_vs_slope_scatter_plot.mean(e,:);
end
plot([0.78 1 1.22],stridelength_variation_detrended_vs_slopes(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
plot([1.78 2 2.22],stridelength_variation_detrended_vs_slopes(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
plot([2.78 3 3.22],stridelength_variation_detrended_vs_slopes(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])

end 
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
ylabel('Step Length Variability -- Detrended vs Slopes')
ylim([0 0.005])

figure 
stridelength_variation_speedtrend_vs_slopes = [mean(variation_steps.decline_050.speedtrend(:)),mean(variation_steps.level_050.speedtrend(:)),mean(variation_steps.incline_050.speedtrend(:));mean(variation_steps.decline_100.speedtrend(:)),mean(variation_steps.level_100.speedtrend(:)),mean(variation_steps.incline_100.speedtrend(:));mean(variation_steps.decline_150.speedtrend(:)),mean(variation_steps.level_150.speedtrend(:)),mean(variation_steps.incline_150.speedtrend(:))];
% bar(stridelength_variation_speedtrend_vs_slopes)
% set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
for i=1:s
stridelength_variation_speedtrend_vs_slope_scatter_plot.mean = [mean(variation_steps.decline_050.speedtrend(i,:),2),mean(variation_steps.level_050.speedtrend(i,:),2),mean(variation_steps.incline_050.speedtrend(i,:),2);mean(variation_steps.decline_100.speedtrend(i,:),2),mean(variation_steps.level_100.speedtrend(i,:),2),mean(variation_steps.incline_100.speedtrend(i,:),2);mean(variation_steps.decline_150.speedtrend(i,:),2),mean(variation_steps.level_150.speedtrend(i,:),2),mean(variation_steps.incline_150.speedtrend(i,:),2)];
plot([0.78 1 1.22],stridelength_variation_speedtrend_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',[0,0,0]+0.8)
plot([1.78 2 2.22],stridelength_variation_speedtrend_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',[0,0,0]+0.8)
plot([2.78 3 3.22],stridelength_variation_speedtrend_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',[0,0,0]+0.8)
for e=1:3
stridelength_variation_speedtrend_vs_slope_scatter_plot.(slopes{e})(i,:) = stridelength_variation_speedtrend_vs_slope_scatter_plot.mean(e,:);
end
plot([0.78 1 1.22],stridelength_variation_speedtrend_vs_slopes(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
plot([1.78 2 2.22],stridelength_variation_speedtrend_vs_slopes(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
plot([2.78 3 3.22],stridelength_variation_speedtrend_vs_slopes(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
end 
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
ylabel('Step Length Variability -- Speedtrend vs Slopes')
ylim([0 0.005])