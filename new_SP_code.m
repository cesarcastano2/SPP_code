close all
clearvars
clc

% subjs = {'SM1' 'SM2' 'SM3' 'SM5' 'SM6' 'SM7' 'SM8' 'SM9' 'SM10' 'SM11'};
subjs = {'SM7'};
% subjs = {'SM10'};
conds = {'level_050' 'level_100' 'level_150' 'incline_050' 'incline_100' 'incline_150' 'decline_050' 'decline_100' 'decline_150'};
% conds = {'level_050' 'level_100' 'level_150'};
theta = [-10 0 10]; 
% projfolder='C:\Users\casta\Documents\Selfpaced modes subjects'; %Laptop
% projfolder='D:\selfpaced modes subjects'; % home desktop 
% projfolder='D:\Cesar\Selfpaced Modes Study'; % lab desktop
% projfolder = pwd;
projfolder='E:\SP Modes'; %for new data

color = jet(length(subjs)); 
color_sub = [0 0 1; 1 0 0; 0 0 0;0 0 1; 1 0 0; 0 0 0;0 0 1; 1 0 0; 0 0 0]; 
% color = [0 0 1; 0 1 0; 1 0 0; 1 0 1; 1 1 0; 0 0 0; 0 1 1; 1 0 1];
pitchval = [0 0 0 10 10 10 -10 -10 -10];
slopes = {'decline' 'level' 'incline'};
for s = 1:length(subjs)
    
    for c = 1:length(conds) 
        %% LOAD DATA FILE
        
        % if you want, you can set the file to load here instead of using the gui
        dflow_file = [projfolder '/' subjs{s} '/' subjs{s} '_' conds{c} '0001.txt'];
        dflow_treadmill_file = [projfolder '/' subjs{s} '/' subjs{s} '_' conds{c} '_treadmill0001.txt'];
        
        if ~exist('dflow_file','var')
            [FileName,PathName,FilterIndex]  = uigetfile('*.txt');
            dflow_file = [PathName FileName];
        end
        
        if ~exist('dflow_treadmill_file','var')
            [FileName,PathName,FilterIndex]  = uigetfile('*treadmill0001.txt');
            dflow_treadmill_file = [PathName FileName];
        end
        
     
        [Frame_df, Time_df, markers_df, forces_df, startidx, stopidx, Total] = import_dflow(dflow_file);
        
        datatreadmill = import_dflow_treadmill(dflow_treadmill_file);

        
        full_data_treadmill = datatreadmill; 
  
       %% Calc Speed from treadmill

       datatreadmill.Time = interpft(datatreadmill.Time(:,:),length(Frame_df));
       datatreadmill.LeftBeltSpeed = interpft(datatreadmill.LeftBeltSpeed(:,:),length(Frame_df));
       datatreadmill.RightBeltSpeed  = interpft(datatreadmill.RightBeltSpeed(:,:),length(Frame_df));
       datatreadmill.Pitch = interpft(datatreadmill.Pitch(:,:),length(Frame_df));
       
       datatreadmill.Time = datatreadmill.Time(startidx:stopidx);
       datatreadmill.LeftBeltSpeed =  datatreadmill.LeftBeltSpeed(startidx:stopidx);
       datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed(startidx:stopidx);
       datatreadmill.Pitch = datatreadmill.Pitch(startidx:stopidx);

       for i=1:length(markers_df.labels)
           markers_df.(markers_df.labels{i}) = markers_df.(markers_df.labels{i})(startidx:stopidx,:);
       end

       forces_df.FP1Cop = forces_df.FP1Cop(startidx:stopidx,:);
       forces_df.FP1For = forces_df.FP1For(startidx:stopidx,:);
       forces_df.FP1Mom = forces_df.FP1Mom(startidx:stopidx,:);
       forces_df.FP2Cop = forces_df.FP2Cop(startidx:stopidx,:);
       forces_df.FP2For = forces_df.FP2For(startidx:stopidx,:);
       forces_df.FP2Mom = forces_df.FP2Mom(startidx:stopidx,:);
        
       Frame_df = Frame_df(startidx:stopidx,:);
       Time_df = Time_df(startidx:stopidx,:);
       
       
             
        %% HARD START / STOP
        Dflow_size = length(Time_df);
        
        if (s==9 && c==5)             %  9/10 
            start_values = find(datatreadmill.LeftBeltSpeed>0.5,1);
        else
            start_values = find(datatreadmill.LeftBeltSpeed>0.1,1);
        end
        start_value_pref = start_values;
            
        Stop_Time_min = floor(Time_df(end));
        Stop_Time = find(Time_df > Stop_Time_min);
        Stop_Time = Stop_Time(1);
       
        Stop_Time = Stop_Time(end);
        
        Start_Time = start_value_pref;
      


        for i=1:length(markers_df.labels)
            markers_df.(markers_df.labels{i}) = markers_df.(markers_df.labels{i})(Start_Time:Stop_Time,:);
        end
        
        forces_df.FP1Cop = forces_df.FP1Cop(Start_Time:Stop_Time,:);
        forces_df.FP1For = forces_df.FP1For(Start_Time:Stop_Time,:);
        forces_df.FP1Mom = forces_df.FP1Mom(Start_Time:Stop_Time,:);
        forces_df.FP2Cop = forces_df.FP2Cop(Start_Time:Stop_Time,:);
        forces_df.FP2For = forces_df.FP2For(Start_Time:Stop_Time,:);
        forces_df.FP2Mom = forces_df.FP2Mom(Start_Time:Stop_Time,:);
        
        Frame_df = Frame_df(Start_Time:Stop_Time,:);
        Time_df = Time_df(Start_Time:Stop_Time,:);
        Time_real = Time_df - Time_df(1);
        
        datatreadmill.Time = datatreadmill.Time(Start_Time:Stop_Time,:);
        datatreadmill.LeftBeltSpeed =  datatreadmill.LeftBeltSpeed(Start_Time:Stop_Time,:);
        datatreadmill.RightBeltSpeed  = datatreadmill.RightBeltSpeed(Start_Time:Stop_Time,:); 
        datatreadmill.Pitch = datatreadmill.Pitch(Start_Time:Stop_Time,:);
          
       
        
%          if c<=3
%          fig1=figure(111);
%          plot(Time_real(1:68400),datatreadmill.LeftBeltSpeed(1:68400),'.-','Color',color_sub(c,:));
%          xlim([0 290])
%          ylim([0 2])
%          gcf = fig1;
%          fig1.PaperUnits = 'inches';
%          fig1.PaperPosition = [0 0 6 3];
%          ylim([0 1.5])
%          hold on
%          end
%          
%          if c==4||c==5||c==6
%          fig2=figure(112);
%          plot(Time_real(1:68400),datatreadmill.LeftBeltSpeed(1:68400),'.-','Color',color_sub(c-3,:));
%          xlim([0 290])
%          ylim([0 2])
%          gcf = fig2;
%          fig2.PaperUnits = 'inches';
%          fig2.PaperPosition = [0 0 6 3];
%          ylim([0 1.5])
%          hold on
%          end
%          
%          if c>=7
%          fig3=figure(113);
%          plot(Time_real(1:68400),datatreadmill.LeftBeltSpeed(1:68400),'.-','Color',color_sub(c-6,:));
%          xlim([0 290])
%          ylim([0 2])
%          gcf = fig3;
%          fig3.PaperUnits = 'inches';
%          fig3.PaperPosition = [0 0 6 3];
%          ylim([0 1.5])
%          hold on
%          end
%         
%       
%     end
% end
        %% GET GAIT EVENTS
        HSrefinePre=10;
        HSrefinePost= 5;
        TOminpeakheight=4 ;
        TOminpeakdistance=40;
        BW=0;
        Time = Time_df;
        markers4GE = {'RHEE' 'LHEE' 'RANK' 'LANK' 'RTOE' 'LTOE'};
%         markers_df.labels(1,:) = {'LASIS', 'RASIS', 'LPSIS', 'RPSIS', 'LKNE', 'LTHI', 'LANK', 'LTIB', 'LTOE', 'LHEE', 'RKNE', 'RTHI', 'RANK', 'RTIB', 'RTOE', 'RHEE'};
        llmarkers= markers_df;
        
        
        [RHS,LTO,LHS,RTO,GE,GEInMiddleDeleted] = GaitEvents_allslopes(Time,llmarkers,markers4GE, HSrefinePre, HSrefinePost, TOminpeakdistance, TOminpeakheight,BW); 
   
        
        RHS = RHS';
        LTO = LTO';
        LHS = LHS';
        RTO = RTO';
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
        for m = 1:length(markers_df.labels)
            eval(['markers_df_c.' markers_df.labels{m} ' = convert_coords2conventional(markers_df. ' markers_df.labels{m} ');']);
        end
        
        for m = 1:length(forces_df.labels)
            eval(['forces_df_c.' forces_df.labels{m} ' = convert_coords2conventional(forces_df. ' forces_df.labels{m} ');']);
        end
        
        %% Filter data, zero-lag fourth order low pass Butterworth filter at 6 Hz
      

        fc = 6;
        fs = 240;
        [b,a] = butter(4,fc/(fs/2));
        full_treadmill_speed.([conds{c} '_filtered'])(:,s)=filtfilt(b,a,datatreadmill.LeftBeltSpeed(1:64000));
        full_treadmill_speed.(conds{c})(:,s) = datatreadmill.LeftBeltSpeed(1:64000);  %R2016b
        
        center_position.(subjs{s}).(conds{c})=(markers_df_c.LASI(:,2)+markers_df_c.RASI(:,2)+markers_df_c.LPSI(:,2)+markers_df_c.RPSI(:,2))/4;

%         datatreadmill.LeftBeltSpeed_plusmarker = datatreadmill.LeftBeltSpeed(1:end-1,1) + diff((markers_df_c.LASI(:,2)+markers_df_c.RASI(:,2)+markers_df_c.LPSI(:,2)+markers_df_c.RPSI(:,2))/4);
%         datatreadmill.RightBeltSpeed_plusmarker = datatreadmill.RightBeltSpeed(1:end-1,1) + diff((markers_df_c.LASI(:,2)+markers_df_c.RASI(:,2)+markers_df_c.LPSI(:,2)+markers_df_c.RPSI(:,2))/4);
%         
%         datatreadmill.RightBeltSpeed_plusmarker = filtfilt(b,a,datatreadmill.RightBeltSpeed_plusmarker);
%         datatreadmill.LeftBeltSpeed_plusmarker = filtfilt(b,a,datatreadmill.LeftBeltSpeed_plusmarker);
        datatreadmill.LeftBeltSpeed_plusmarker = datatreadmill.LeftBeltSpeed(1:end-1,1) ;
        datatreadmill.RightBeltSpeed_plusmarker = datatreadmill.RightBeltSpeed(1:end-1,1) ;

%       markers_df_c_filtered = [];
        for i=1:length(markers_df.labels)
            for e=1:3
                markers_df_c.(markers_df.labels{i})(:,e) = filtfilt(b,a,markers_df_c.(markers_df.labels{i})(:,e));
            end
        end

        allmarkers_c_filtered.(subjs{s}).(conds{c}) = markers_df_c;
        allmarkers_c.(subjs{s}).(conds{c}) = markers_df_c;

        % mark_plus_treadspeed=diff(markers_df_filtered.RASI(:,2))+datatreadmill.RightBeltSpeed(1:(end-1),1);

        %% step length / Step Time / Stride Time
        stride_speed = [];
        stride_time = [];
        steplength_speed = [];
        steplength_time = [];
        stridelength_value = [];
        stridewidth_value = [];
        vertical_GRF = [];
       
        
        
        for i = 1:length(GEgood(:,1))
            vertical_GRF(i,1) = nanmean(forces_df_c.FP2For((GEgood(i,1):GEgood(i,3)),3));
            vertical_GRF(i,2) = nanmean(forces_df_c.FP1For((GEgood(i,3):GEgood(i,5)),3));
        end  

        for i = 1:length(GEgood(:,1))
            steplength_speed(i,1) = nanmean(datatreadmill.RightBeltSpeed_plusmarker(GEgood(i,1):GEgood(i,3)));
            steplength_speed(i,2) = nanmean(datatreadmill.LeftBeltSpeed_plusmarker(GEgood(i,3):GEgood(i,5)));
        end  

        for i = 1:length(GEgood(:,1))
            steplength_time(i,1) = range(Time_df(GEgood(i,1):GEgood(i,3)));
            steplength_time(i,2) = range(Time_df(GEgood(i,3):GEgood(i,5)));
        end
        
        for i = 1:length(GEgood(:,1))
            stride_time(i,1) = range(Time_df(GEgood(i,1):GEgood(i,5)));
            stride_time(i,2) = range(Time_df(GEgood(i,3):leftheel_GE(i,1)));
        end
        
      
        for i = 1:length(GEgood(:,1))
            stride_speed(i,1) = nanmean(datatreadmill.RightBeltSpeed_plusmarker(GEgood(i,1):GEgood(i,5)));
            stride_speed(i,2) = nanmean(datatreadmill.LeftBeltSpeed_plusmarker(GEgood(i,3):leftheel_GE(i,1)));    
        end
        
        stridelength_value(:,1) = stridelength(markers_df_c.RHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,5), GEgood(:,1),stride_speed(:,1), stride_time(:,1),pitchval(c));
        stridelength_value(:,2) = stridelength(markers_df_c.LHEE(:,2), markers_df_c.LHEE(:,2), leftheel_GE(:,1), GEgood(:,3),stride_speed(:,2), stride_time(:,2),pitchval(c));
        
        stridewidth_value(:,1) = stridewidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,1), GEgood(:,3));
        stridewidth_value(:,2) = abs(stridewidth(markers_df_c.LHEE(:,1),markers_df_c.RHEE(:,1),  GEgood(:,3), GEgood(:,5)));
        %%
%         h(c)=figure;
%         GEgood = gaiteventCheck3(RHS, LTO, LHS, RTO, forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3));

        figure
        
        plot(forces_df_c.FP2For(:,3),'r');
        hold on
        plot(forces_df_c.FP1For(:,3),'b');
        %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k', 'linewidth', 2)
        %         plot(min(forces_df_c.FP2For(:,3), forces_df_c.FP1For(:,3)), 'k')
        plot(GEgood(:,1), zeros(size(GEgood(:,1))), 'rx', GEgood(:,2), zeros(size(GEgood(:,1))), 'bo', GEgood(:,3), zeros(size(GEgood(:,1))), 'bx', GEgood(:,4), zeros(size(GEgood(:,1))), 'ro', GEgood(:,5), zeros(size(GEgood(:,1))), 'rx')
        ylimits = ylim(gca);
        title([subjs{s} ', ' conds{c}])
        %             axis([0 RHS(10) ylimits(1) ylimits(2)]);
        
        
%         savefig(h,'Incline_gait_ACC.fig')
%         close(gcf)
        %% SAVE GE VARIABLES
        
%         filename = [subjs(s) conds(c)];
%         filename = strjoin(filename,{'_'});
%         filename =  strcat(filename,'.mat');
%         save(filename,'GE')
        
        %% Step width & Step length
        
        clear sl sw badidx1 badidx2 allbadsl
        
        sl=[]; sw=[]; sp=[];  st=[];
        
        sl(:,1) = steplength(markers_df_c.LHEE(:,2), markers_df_c.RHEE(:,2), GEgood(:,3), GEgood(:,1),steplength_speed(:,1), steplength_time(:,1),pitchval(c));
        sw(:,1) = stepwidth(markers_df_c.RHEE(:,1), markers_df_c.LHEE(:,1), GEgood(:,3));
        
        sl(:,2) = steplength(markers_df_c.RHEE(:,2), markers_df_c.LHEE(:,2), GEgood(:,5), GEgood(:,3),steplength_speed(:,2), steplength_time(:,2),pitchval(c));
        sw(:,2) = stepwidth(markers_df_c.LHEE(:,1), markers_df_c.RHEE(:,1), GEgood(:,5));
        goodidx1 = find(sl(:,1) > 0.9*mean(sl(:,1)));
        goodidx2 = find(sl(:,2) > 0.9*mean(sl(:,2)));
        goodsl = intersect(goodidx1, goodidx2); 
        tempsl = sl(goodsl,:);
        tempsw = sw(goodsl,:);
        stepcount = length(sl);
%         sl = tempsl;
%         sw = tempsw;
        sp = datatreadmill.LeftBeltSpeed_plusmarker(:);
%         sp = avespeed(39:end)';
        for sn = 1:length(GEgood)
            st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
        end
%         for sn = 41:length(GEgood)
%             st(sn,1) = Time_df(GEgood(sn,5))-Time_df(GEgood(sn,1));
%         end
%         
        GRP.sp(s,c) = mean(sp);
        GRP.st(s,c) = mean(st);
        GRP.sl(s,:,c) = mean(sl);
        GRP.sw(s,:,c) = mean(sw);
        
        GRP.spvar(s,c) = std(sp);
        GRP.stvar(s,c) = std(st);
        GRP.slvar(s,:,c) = std(sl);
        GRP.swvar(s,:,c) = std(sw);
        %% GRF
        
        numstrides = size(GEgood,1);
        npts = 250;
        GEn = 1; %RHS
        
        [time_n, GRF_n, data_nft] = normalize_gaitcycle2(Time_df, forces_df_c.FP2For, GEgood, npts);
        
        temp = nanmean(GRF_n,3);
        GRP.GRF_ML_R(:,s,c) = temp(:,1);
        GRP.GRF_AP_R(:,s,c) = temp(:,2);
        GRP.GRF_v_R(:,s,c) = temp(:,3);            
%% PLOTS VS STEPS
 
sl_all = [];
for i = 1:size(sl,1)
    sl_all = [sl_all sl(i,:)];
end

sw_all = [];
for i = 1:size(sw,1)
    sw_all = [sw_all sw(i,:)];
end

% steplength_time = steplength_time';
steplength_time_all = [];
for i = 1:size(steplength_time,1)
    steplength_time_all = [steplength_time_all steplength_time(i,:)];
end

vertical_GRF_all =[];
for i = 1:size(vertical_GRF,1)
    vertical_GRF_all = [vertical_GRF_all vertical_GRF(i,:)];
end
% steplength_speed = steplength_speed';
steplength_speed_all = [];
for i = 1:size(steplength_speed,1)
    steplength_speed_all = [steplength_speed_all steplength_speed(i,:)];
end

stride_speed_all = [];
for i = 1:size(stride_speed,1)
    stride_speed_all = [stride_speed_all stride_speed(i,:)];
end

stride_time_all = [];
for i = 1:size(stride_time,1)
    stride_time_all = [stride_time_all stride_time(i,:)];
end

stridelength_all = [];
for i = 1:size(stridelength_value,1)
    stridelength_all = [stridelength_all stridelength_value(i,:)];
end

stridewidth_all = [];
for i = 1:size(stridewidth_value,1)
    stridewidth_all = [stridewidth_all stridewidth_value(i,:)];
end  
%% STD Method 
% [STV_OP1_value,STV_OP1_median_value,STV_OP1_steadypoint_value_all,STV_OP1_steadypoint_value] = stdmethod(stride_speed_all(1:375));
% 
% if c == 1
%   STV_OP1.level.steadystatevalue(:,s) = STV_OP1_steadypoint_value;
%   STV_OP1.level.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/269;
%   STV_OP1.level.medianvalue(:,s) = STV_OP1_median_value;
%   STV_OP1.level.STDvalues(:,s) = STV_OP1_value;
% end
% if c == 2      
%   STV_OP1.incline.steadystatevalue(:,s) = STV_OP1_steadypoint_value;
%   STV_OP1.incline.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/269;
%   STV_OP1.incline.medianvalue(:,s) = STV_OP1_median_value;
%   STV_OP1.incline.STDvalues(:,s) = STV_OP1_value;
% end
% if c == 3
%   STV_OP1.decline.steadystatevalue(:,s) = STV_OP1_steadypoint_value; 
%   STV_OP1.decline.steadystatevalue_all(:,s) = length(STV_OP1_steadypoint_value_all)/269;
%   STV_OP1.decline.medianvalue(:,s) = STV_OP1_median_value;
%   STV_OP1.decline.STDvalues(:,s) = STV_OP1_value;
% end
% STV_OP1_steadypoint_value = [];
% STV_OP1_value = [];
% STV_OP1_median_value = [];
% STV_OP1_steadypoint_value_all = [];

%% Speed plots

% figure(110)
% subplot(3,3,c)
% title([subjs{s} ' ' conds{c}])
% plot(1:length(stride_speed_all(1:269)),stride_speed_all(1:269),'.-','Color',color(s,:)); 
% xlim([0 269])
% hold on 
% % 


% figure(111)
% subplot(3,3,c)
% title([subjs{s} ' ' conds{c}])
% plot(1:length(datatreadmill.LeftBeltSpeed(:)),datatreadmill.LeftBeltSpeed(:),'.-','Color',color(s,:)); 
% xlim([0 65000])
% ylim([0 2])
% hold on 
% 
%% Power Spectral Density 

PSD_window=(find(Time_real>60,1):find(Time_real>60,1)+43203)';
Speed_PSD_without_mean=datatreadmill.LeftBeltSpeed_plusmarker(PSD_window);  
Speed_PSD_mean = Speed_PSD_without_mean - mean(Speed_PSD_without_mean);
[Power_PSD, Frequency_PSD] = powerspectrumplot(Speed_PSD_mean);


 PowerSpectrumValues.([conds{c} '_PSD'])(:,s) = Power_PSD;
 PowerSpectrumValues.([conds{c} '_PowerPSD'])(:,s) = Power_PSD./sum(Power_PSD);

%% Stride speeds, lengths, frequency for all subjects 
%GRF 
GRF_all_subjs.(conds{c}).(subjs{s})= vertical_GRF_all;
%Step Speed 
speed_steps.(conds{c}).(subjs{s}) = steplength_speed_all; 
% Step Length 
length_steps.(conds{c}).(subjs{s}) = sl_all; 
% Step Width
width_steps.(conds{c}).(subjs{s}) = sw_all; 
% Step Time 
frequency_steps.(conds{c}).(subjs{s}) = 1./steplength_time_all;
% % Subject Speed 
% speed_SLandST.(conds{c})(s,:) = sl_all(1:375)./steplength_time_all(1:375);
% % Stride Speed 
% speed_strides.(conds{c})(s,:) = stride_speed_all(1:275);
% % Stride Length  
% length_strides.(conds{c})(s,:) = stridelength_all(1:275);     
% % Stride Width  
% width_strides.(conds{c})(s,:) = stridewidth_all(1:275);
% % Stride Time  
% time_strides.(conds{c})(s,:) = stride_time_all(1:275);    
%% for fixed speeds
% p_fixed.(subjs{s})(c,:) = P;
NOG = (length(steplength_speed_all)) - 325;

if c<=3
length_steps_forfixedthing.(subjs{s}).all_level_conds(c,:) = length_steps.(conds{c}).(subjs{s})(NOG:end);
end
if c==4||c==5||c==6
length_steps_forfixedthing.(subjs{s}).all_incline_conds(c-3,:) = length_steps.(conds{c}).(subjs{s})(NOG:end);
end
if c>=7
length_steps_forfixedthing.(subjs{s}).all_decline_conds(c-6,:) = length_steps.(conds{c}).(subjs{s})(NOG:end);
end

if c<=3
speed_steps_forfixedthing.(subjs{s}).all_level_conds(c,:) = speed_steps.(conds{c}).(subjs{s})(NOG:end);
end
if c==4||c==5||c==6
speed_steps_forfixedthing.(subjs{s}).all_incline_conds(c-3,:) = speed_steps.(conds{c}).(subjs{s})(NOG:end);
end
if c>=7
speed_steps_forfixedthing.(subjs{s}).all_decline_conds(c-6,:) = speed_steps.(conds{c}).(subjs{s})(NOG:end);
end


% [P,fh] = fitsinglemodelprocess_sl(length_steps_forfixedthing.(subjs{s}).all_decline_conds(c,:),speed_steps_forfixedthing.(subjs{s}).all_decline_conds(c,:));           
% p.(conds{c})(s,:) = P; 
% p_fixed.(subjs{s})(c,:) = P;
%% FINDING VARIABILITY FOR STRIDE LENGTH 

% [P,fh] = fitsinglemodelprocess_sl(length_steps.(conds{c})(s,75:375),speed_steps.(conds{c})(s,75:375));           
% p.(conds{c})(s,:) = P;   
% 
% variation_steps_width.(conds{c}).slminusfit(s,:) = var(length_steps.(conds{c})(s,75:375) - fh(speed_steps.(conds{c})(s,75:375),P));
% variation_steps.(conds{c}).speedtrend(s,:) = var(fh(speed_steps.(conds{c})(s,75:375),P));
% variation_steps.(conds{c}).totalvar(s,:) = var(length_steps.(conds{c})(s,75:375));
% variation_steps.(conds{c}).totalstd(s,:) = std(length_steps.(conds{c})(s,75:375));
% variation_steps.(conds{c}).totalcov(s,:) = cov(length_steps.(conds{c})(s,75:375));
% variation_steps.(conds{c}).speedvar(s,:) = std(speed_steps.(conds{c})(s,75:375));

NOG = (length(steplength_speed_all)) - 325;


[P,fh] = fitsinglemodelprocess_sl(length_steps.(conds{c}).(subjs{s})(NOG:end),speed_steps.(conds{c}).(subjs{s})(NOG:end));           
p.(conds{c})(s,:) = P; 

fitplot_fitted_steplength.(conds{c}).(subjs{s}) = fh(speed_steps.(conds{c}).(subjs{s})(NOG:end),P);
fitplot_stepfitted_minus_actualstep.(conds{c}).(subjs{s}) = (length_steps.(conds{c}).(subjs{s})(NOG:end) - fh(speed_steps.(conds{c}).(subjs{s})(NOG:end),P));
fitplot_speed.(conds{c}).(subjs{s}) = speed_steps.(conds{c}).(subjs{s})(NOG:end);
fitplot_actual_steplength.(conds{c}).(subjs{s}) = length_steps.(conds{c}).(subjs{s})(NOG:end);


variation_steps.(conds{c}).slminusfit(s,:) = var(length_steps.(conds{c}).(subjs{s})(NOG:end) - fh(speed_steps.(conds{c}).(subjs{s})(NOG:end),P));
variation_steps.(conds{c}).speedtrend(s,:) = var(fh(speed_steps.(conds{c}).(subjs{s})(NOG:end),P));
variation_steps.(conds{c}).totalvar(s,:) = var(length_steps.(conds{c}).(subjs{s})(NOG:end));
variation_steps.(conds{c}).totalstd(s,:) = std(length_steps.(conds{c}).(subjs{s})(NOG:end));
variation_steps.(conds{c}).totalcov(s,:) = cov(length_steps.(conds{c}).(subjs{s})(NOG:end));
variation_steps.(conds{c}).speedvar(s,:) = var(speed_steps.(conds{c}).(subjs{s})(NOG:end));
variation_steps.(conds{c}).speedmean(s,:) = mean(speed_steps.(conds{c}).(subjs{s})(NOG:end));
length_steps.(conds{c}).steplength_mean(s) = mean(length_steps.(conds{c}).(subjs{s})(NOG:end)); 
length_steps.(conds{c}).steplength_var(s) = var(length_steps.(conds{c}).(subjs{s})(NOG:end)); 
width_steps.(conds{c}).stepwidth_mean(s) = mean(width_steps.(conds{c}).(subjs{s})(NOG:end)); 
width_steps.(conds{c}).stepwidth_var(s) = var(width_steps.(conds{c}).(subjs{s})(NOG:end)); 
frequency_steps.(conds{c}).stepfreq_mean(s) = mean(frequency_steps.(conds{c}).(subjs{s})(NOG:end));
frequency_steps.(conds{c}).stepfreq_var(s) = var(frequency_steps.(conds{c}).(subjs{s})(NOG:end));
walking_speed.(conds{c}).walkingspeed_mean(s) = mean(speed_steps.(conds{c}).(subjs{s})(NOG:end));
walking_speed.(conds{c}).walkingspeed_var(s) = var(speed_steps.(conds{c}).(subjs{s})(NOG:end));
%% FINDING VARIABILITY FOR STRIDE WIDTH
% 
% [P_width,fh_width] = fitsinglemodelprocess_sw(width_steps.(conds{c})(s,75:375),speed_steps.(conds{c})(s,75:375));
% p_width.(conds{c})(s,:) = P_width;
%  
% variation_steps_width.(conds{c}).slminusfit(s,:) = var(width_steps.(conds{c})(s,75:375) - fh_width(speed_steps.(conds{c})(s,75:375),P_width));
% variation_steps_width.(conds{c}).speedtrend(s,:) = var(fh_width(speed_steps.(conds{c})(s,75:375),P_width));
% variation_steps_width.(conds{c}).totalvar(s,:) = var(width_steps.(conds{c})(s,75:375));


[P_width,fh_width] = fitsinglemodelprocess_sw(width_steps.(conds{c}).(subjs{s})(NOG:end),speed_steps.(conds{c}).(subjs{s})(NOG:end));
p_width.(conds{c})(s,:) = P_width;
 
variation_steps_width.(conds{c}).slminusfit(s,:) = var(width_steps.(conds{c}).(subjs{s})(NOG:end) - fh_width(speed_steps.(conds{c}).(subjs{s})(NOG:end),P_width));
variation_steps_width.(conds{c}).speedtrend(s,:) = var(fh_width(speed_steps.(conds{c}).(subjs{s})(NOG:end),P_width));
variation_steps_width.(conds{c}).totalvar(s,:) = var(width_steps.(conds{c}).(subjs{s})(NOG:end));
%% TREADMILL POSITION
% center_position.(subjs{s}).(conds{c})=(markers_df_c.LASI(:,2)+markers_df_c.RASI(:,2)+markers_df_c.LPSI(:,2)+markers_df_c.RPSI(:,2))/4;
    end 
end
%% fixed stuff 
% p_fixed.(subjs{s})(c,:) = P;
% NOG = (length(steplength_speed_all)) - 325;
% 
% for m=1:s
%         length_steps_forfixedthing.(subjs{m}).all_level_conds_combined(1,:)=[length_steps_forfixedthing.(subjs{m}).all_level_conds(1,:) length_steps_forfixedthing.(subjs{m}).all_level_conds(2,:) length_steps_forfixedthing.(subjs{m}).all_level_conds(3,:)];  
%         length_steps_forfixedthing.(subjs{m}).all_incline_conds_combined(1,:)=[length_steps_forfixedthing.(subjs{m}).all_incline_conds(1,:) length_steps_forfixedthing.(subjs{m}).all_incline_conds(2,:) length_steps_forfixedthing.(subjs{m}).all_incline_conds(3,:)];
%         length_steps_forfixedthing.(subjs{m}).all_decline_conds_combined(1,:)=[length_steps_forfixedthing.(subjs{m}).all_decline_conds(1,:) length_steps_forfixedthing.(subjs{m}).all_decline_conds(2,:) length_steps_forfixedthing.(subjs{m}).all_decline_conds(3,:)];
%         speed_steps_forfixedthing.(subjs{m}).all_level_conds_combined(1,:)=[speed_steps_forfixedthing.(subjs{m}).all_level_conds(1,:) speed_steps_forfixedthing.(subjs{m}).all_level_conds(2,:) speed_steps_forfixedthing.(subjs{m}).all_level_conds(3,:)];  
%         speed_steps_forfixedthing.(subjs{m}).all_incline_conds_combined(1,:)=[speed_steps_forfixedthing.(subjs{m}).all_incline_conds(1,:) speed_steps_forfixedthing.(subjs{m}).all_incline_conds(2,:) speed_steps_forfixedthing.(subjs{m}).all_incline_conds(3,:)];
%         speed_steps_forfixedthing.(subjs{m}).all_decline_conds_combined(1,:)=[speed_steps_forfixedthing.(subjs{m}).all_decline_conds(1,:) speed_steps_forfixedthing.(subjs{m}).all_decline_conds(2,:) speed_steps_forfixedthing.(subjs{m}).all_decline_conds(3,:)];
% 
% end
% 
% for m=1:s
% [P_level,fh] = fitsinglemodelprocess_sl(length_steps_forfixedthing.(subjs{m}).all_level_conds_combined(1,:), speed_steps_forfixedthing.(subjs{m}).all_level_conds_combined(1,:));           
% [P_incline,fh] = fitsinglemodelprocess_sl(length_steps_forfixedthing.(subjs{m}).all_incline_conds_combined(1,:), speed_steps_forfixedthing.(subjs{m}).all_incline_conds_combined(1,:));           
% [P_decline,fh] = fitsinglemodelprocess_sl(length_steps_forfixedthing.(subjs{m}).all_decline_conds_combined(1,:), speed_steps_forfixedthing.(subjs{m}).all_decline_conds_combined(1,:));           
% p_fixed.level.(subjs{m})=P_level;
% p_fixed.incline.(subjs{m})=P_incline;
% p_fixed.decline.(subjs{m})=P_decline;
% end

%% VARIANCE PLOTS
% detrended 
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_detrended(1,i-6)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (1<=i&&i<=3)
        error_length_detrended(2,i)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (4<=i&&i<=6)
        error_length_detrended(3,i-3)=std(variation_steps.(conds{i}).slminusfit(:));
    end 
end
figure 
subplot(1,3,2)
stridelength_variation_detrended_bar_plot = [mean(variation_steps.decline_050.slminusfit(:)),mean(variation_steps.decline_100.slminusfit(:)),mean(variation_steps.decline_150.slminusfit(:));mean(variation_steps.level_050.slminusfit(:)),mean(variation_steps.level_100.slminusfit(:)),mean(variation_steps.level_150.slminusfit(:));mean(variation_steps.incline_050.slminusfit(:)),mean(variation_steps.incline_100.slminusfit(:)),mean(variation_steps.incline_150.slminusfit(:))];
bar(stridelength_variation_detrended_bar_plot);
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_detrended_bar_plot,error_length_detrended,'.');
for i=1:s
stridelength_variation_detrended_scatter_plot.mean = [mean(variation_steps.decline_050.slminusfit(i,:),2),mean(variation_steps.decline_100.slminusfit(i,:),2),mean(variation_steps.decline_150.slminusfit(i,:),2);mean(variation_steps.level_050.slminusfit(i,:),2),mean(variation_steps.level_100.slminusfit(i,:),2),mean(variation_steps.level_150.slminusfit(i,:),2);mean(variation_steps.incline_050.slminusfit(i,:),2),mean(variation_steps.incline_100.slminusfit(i,:),2),mean(variation_steps.incline_150.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_detrended_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_detrended_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_detrended_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_detrended_scatter_plot.(slopes{e})(i,:) = stridelength_variation_detrended_scatter_plot.mean(e,:);
end
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
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_detrended(1,i-6)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (1<=i&&i<=3)
        error_length_detrended(2,i)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (4<=i&&i<=6)
        error_length_detrended(3,i-3)=std(variation_steps.(conds{i}).slminusfit(:));
    end 
end
stridelength_stacked=[];
for i=[1 4 7 2 5 8 3 6 9]
e=size(stridelength_stacked,1)+1;
stridelength_stacked(e,:)=[stridelength_variation_detrended_bar_plot(i),stridelength_variation_speedtrend_bar_plot(i)];
stridelength_stacked(e+1,:)=[stridelength_variation_total_bar_plot(i),0]; 
end 
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
%% Stride Width Variability 
% Detrended 
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_width_detrended(1,i-6)=std(variation_steps_width.(conds{i}).slminusfit(:));
    elseif (1<=i&&i<=3)
        error_width_detrended(2,i)=std(variation_steps_width.(conds{i}).slminusfit(:));
    elseif (4<=i&&i<=6)
        error_width_detrended(3,i-3)=std(variation_steps_width.(conds{i}).slminusfit(:));
    end 
end
figure 
subplot(1,3,2)
stridewidth_variation_detrended_bar_plot = [mean(variation_steps_width.decline_050.slminusfit(:)),mean(variation_steps_width.decline_100.slminusfit(:)),mean(variation_steps_width.decline_150.slminusfit(:));mean(variation_steps_width.level_050.slminusfit(:)),mean(variation_steps_width.level_100.slminusfit(:)),mean(variation_steps_width.level_150.slminusfit(:));mean(variation_steps_width.incline_050.slminusfit(:)),mean(variation_steps_width.incline_100.slminusfit(:)),mean(variation_steps_width.incline_150.slminusfit(:))];
bar(stridewidth_variation_detrended_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
for i=1:s
stridewidth_variation_detrended_scatter_plot.mean = [mean(variation_steps_width.decline_050.slminusfit(i,:),2),mean(variation_steps_width.decline_100.slminusfit(i,:),2),mean(variation_steps_width.decline_150.slminusfit(i,:),2);mean(variation_steps_width.level_050.slminusfit(i,:),2),mean(variation_steps_width.level_100.slminusfit(i,:),2),mean(variation_steps_width.level_150.slminusfit(i,:),2);mean(variation_steps_width.incline_050.slminusfit(i,:),2),mean(variation_steps_width.incline_100.slminusfit(i,:),2),mean(variation_steps_width.incline_150.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_detrended_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_detrended_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_detrended_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_detrended_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_detrended_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Detrended')
ylim([0 0.0006])
% Total Variability 
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_width_totalvar(1,i-6)=std(variation_steps_width.(conds{i}).totalvar(:));
    elseif (1<=i&&i<=3)
        error_width_totalvar(2,i)=std(variation_steps_width.(conds{i}).totalvar(:));
    elseif (4<=i&&i<=6)
        error_width_totalvar(3,i-3)=std(variation_steps_width.(conds{i}).totalvar(:));
    end 
end
subplot(1,3,3) 
stridewidth_variation_total_bar_plot = [mean(variation_steps_width.decline_050.totalvar(:)),mean(variation_steps_width.decline_100.totalvar(:)),mean(variation_steps_width.decline_150.totalvar(:));mean(variation_steps_width.level_050.totalvar(:)),mean(variation_steps_width.level_100.totalvar(:)),mean(variation_steps_width.level_150.totalvar(:));mean(variation_steps_width.incline_050.totalvar(:)),mean(variation_steps_width.incline_100.totalvar(:)),mean(variation_steps_width.incline_150.totalvar(:))];
bar(stridewidth_variation_total_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridewidth_variation_total_bar_plot,error_width_totalvar,'.');
for i=1:s
stridewidth_variation_total_scatter_plot.mean = [mean(variation_steps_width.decline_050.totalvar(i,:),2),mean(variation_steps_width.decline_100.totalvar(i,:),2),mean(variation_steps_width.decline_150.totalvar(i,:),2);mean(variation_steps_width.level_050.totalvar(i,:),2),mean(variation_steps_width.level_100.totalvar(i,:),2),mean(variation_steps_width.level_150.totalvar(i,:),2);mean(variation_steps_width.incline_050.totalvar(i,:),2),mean(variation_steps_width.incline_100.totalvar(i,:),2),mean(variation_steps_width.incline_150.totalvar(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_total_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_total_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Total')
ylim([0 0.0006])
% Speedtrend
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_width_speedtrend(1,i-6)=std(variation_steps_width.(conds{i}).speedtrend(:));
    elseif (1<=i&&i<=3)
        error_width_speedtrend(2,i)=std(variation_steps_width.(conds{i}).speedtrend(:));
    elseif (4<=i&&i<=6)
        error_width_speedtrend(3,i-3)=std(variation_steps_width.(conds{i}).speedtrend(:));
    end 
end
subplot(1,3,1)
stridewidth_variation_speedtrend_bar_plot = [mean(variation_steps_width.decline_050.speedtrend(:)),mean(variation_steps_width.decline_100.speedtrend(:)),mean(variation_steps_width.decline_150.speedtrend(:));mean(variation_steps_width.level_050.speedtrend(:)),mean(variation_steps_width.level_100.speedtrend(:)),mean(variation_steps_width.level_150.speedtrend(:));mean(variation_steps_width.incline_050.speedtrend(:)),mean(variation_steps_width.incline_100.speedtrend(:)),mean(variation_steps_width.incline_150.speedtrend(:))];
bar(stridewidth_variation_speedtrend_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridewidth_variation_speedtrend_bar_plot,error_width_speedtrend,'.');
for i=1:s
stridewidth_variation_speedtrend_scatter_plot.mean = [mean(variation_steps_width.decline_050.speedtrend(i,:),2),mean(variation_steps_width.decline_100.speedtrend(i,:),2),mean(variation_steps_width.decline_150.speedtrend(i,:),2);mean(variation_steps_width.level_050.speedtrend(i,:),2),mean(variation_steps_width.level_100.speedtrend(i,:),2),mean(variation_steps_width.level_150.speedtrend(i,:),2);mean(variation_steps_width.incline_050.speedtrend(i,:),2),mean(variation_steps_width.incline_100.speedtrend(i,:),2),mean(variation_steps_width.incline_150.speedtrend(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_speedtrend_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_speedtrend_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_speedtrend_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_speedtrend_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_speedtrend_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Speedtrend')
ylim([0 0.0006])
% Detrended vs Slope
figure 
stridewidth_variation_detrended_vs_slopes = [mean(variation_steps_width.decline_050.slminusfit(:)),mean(variation_steps_width.level_050.slminusfit(:)),mean(variation_steps_width.incline_050.slminusfit(:));mean(variation_steps_width.decline_100.slminusfit(:)),mean(variation_steps_width.level_100.slminusfit(:)),mean(variation_steps_width.incline_100.slminusfit(:));mean(variation_steps_width.decline_150.slminusfit(:)),mean(variation_steps_width.level_150.slminusfit(:)),mean(variation_steps_width.incline_150.slminusfit(:))];
bar(stridewidth_variation_detrended_vs_slopes)
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
for i=1:s
stridewidth_variation_detrended_vs_slope_scatter_plot.mean = [mean(variation_steps_width.decline_050.slminusfit(i,:),2),mean(variation_steps_width.level_050.slminusfit(i,:),2),mean(variation_steps_width.incline_050.slminusfit(i,:),2);mean(variation_steps_width.decline_100.slminusfit(i,:),2),mean(variation_steps_width.level_100.slminusfit(i,:),2),mean(variation_steps_width.incline_100.slminusfit(i,:),2);mean(variation_steps_width.decline_150.slminusfit(i,:),2),mean(variation_steps_width.level_150.slminusfit(i,:),2),mean(variation_steps_width.incline_150.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_detrended_vs_slope_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_detrended_vs_slope_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Detrended vs Slopes')
% ylim([0 0.02])
%% Position 
% figure 
% subplot(1,3,1)
% for e=1
%     for i=1:length(subjs)
% plot(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2)), hold on 
%     end
% end  
% subplot(1,3,2)
% for e=2
%     for i=1:length(subjs)
% plot(allmarkers_c.(subjs{i}).(conds{i}).RASI(10000:55000,2)), hold on 
%     end
% end 
% subplot(1,3,3)
% for e=3
%     for i=1:length(subjs)
% plot(allmarkers_c.(subjs{i}).(conds{i}).RASI(10000:55000,2)), hold on 
%     end
% end 

% for e=1:length(conds)
%     for i=1:length(subjs)
%         RASI_position.var.(conds{e})(i,1) = var(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2));
%         RASI_position.std.(conds{e})(i,1) = std(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2));
%         RASI_position.range.(conds{e})(i,1) = range(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2));
%         RASI_velocity.var.(conds{e})(i,1) = var(diff(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2)));
%         RASI_velocity.std.(conds{e})(i,1) = std(diff(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2)));
%     end
%     RASI_position.var.mean(e,1) = mean(RASI_position.var.(conds{e})(:,1));
%     RASI_position.std.mean(e,1) = mean(RASI_position.std.(conds{e})(:,1));
%     RASI_position.range.mean(e,1) = mean(RASI_position.range.(conds{e})(:,1));
%     RASI_velocity.var.mean(e,1) = mean(RASI_velocity.var.(conds{e})(:,1));
%     RASI_velocity.std.mean(e,1) = mean(RASI_velocity.std.(conds{e})(:,1));
% end
% 
% error_position(1,1:3) = RASI_position.std.mean(7:9);
% error_position(2,1:3) = RASI_position.std.mean(1:3);
% error_position(3,1:3) = RASI_position.std.mean(3:5);
% 
% barerror_position(1,1:3) = RASI_position.range.mean(7:9);
% barerror_position(2,1:3) = RASI_position.range.mean(1:3);
% barerror_position(3,1:3) = RASI_position.range.mean(4:6);
% 
% for e=1:length(conds)
%     for i=1:length(subjs)
%         RASI_position.velocity.(conds{e})(i,:) = diff(allmarkers_c.(subjs{i}).(conds{e}).RASI(10000:55000,2));
%     end
% end

% figure 
% bar([RASI_position.range.mean(7),RASI_position.range.mean(8),RASI_position.range.mean(9);RASI_position.range.mean(1),RASI_position.range.mean(2),RASI_position.range.mean(3);RASI_position.range.mean(4),RASI_position.range.mean(5),RASI_position.range.mean(6)])
% hold on
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],barerror_position,error_position,'.');
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% ylabel(['Range of Position (m)'])
%% fun plots 
% x=[(length_steps.level_050.SM1) ; (length_steps.level_050.SM2)]
% y
% table_steplength=table([repmat(subjs(1),1,300),repmat(subjs(2),1,300)]',[length_steps.level_050.SM1(1:300),length_steps.level_050.SM2(1:300)]',[(1:300),(1:300)]','VariableNames',{'subject','step_length','step_number'});
% figure 
% gscatter(table_steplength.step_number,table_steplength.step_length,table_steplength.subject)
% xlabel("step #")
% ylabel("Step length")
%% NEW FIGURE 1 (SPEED)
% fig_2=figure
% subplot(2,4,1) 
% for i=[7 8 9 1 2 3 4 5 6]
%     if (7<=i&&i<=9)
%         err_speedmean(1,i-6)=std(variation_steps.(conds{i}).speedmean(:));
%     elseif (1<=i&&i<=3)
%         err_speedmean(2,i)=std(variation_steps.(conds{i}).speedmean(:));
%     elseif (4<=i&&i<=6)
%         err_speedmean(3,i-3)=std(variation_steps.(conds{i}).speedmean(:));
%     end 
% end
% stride_variation_speedmean_bar_plot = [mean(variation_steps.decline_050.speedmean(:)),mean(variation_steps.decline_100.speedmean(:)),mean(variation_steps.decline_150.speedmean(:));mean(variation_steps.level_050.speedmean(:)),mean(variation_steps.level_100.speedmean(:)),mean(variation_steps.level_150.speedmean(:));mean(variation_steps.incline_050.speedmean(:)),mean(variation_steps.incline_100.speedmean(:)),mean(variation_steps.incline_150.speedmean(:))];
% % barplots.speedmean=bar(stride_variation_speedmean_bar_plot);
% % set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% % hold on 
% % errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stride_variation_speedmean_bar_plot,err_speedmean,'.');
% for i=1:s
% stride_variation_speedmean_scatter_plot.mean = [mean(variation_steps.decline_050.speedmean(i,:),2),mean(variation_steps.decline_100.speedmean(i,:),2),mean(variation_steps.decline_150.speedmean(i,:),2);mean(variation_steps.level_050.speedmean(i,:),2),mean(variation_steps.level_100.speedmean(i,:),2),mean(variation_steps.level_150.speedmean(i,:),2);mean(variation_steps.incline_050.speedmean(i,:),2),mean(variation_steps.incline_100.speedmean(i,:),2),mean(variation_steps.incline_150.speedmean(i,:),2)];
% % plot([0.78 1 1.22],stride_variation_speedmean_scatter_plot.mean(1,:),'Marker','.','Color',[0,0,0]+0.8); hold on
% % plot([1.78 2 2.22],stride_variation_speedmean_scatter_plot.mean(2,:),'Marker','.','Color',[0,0,0]+0.8)
% % plot([2.78 3 3.22],stride_variation_speedmean_scatter_plot.mean(3,:),'Marker','.','Color',[0,0,0]+0.8)
% plot([0.78 1 1.22],stride_variation_speedmean_scatter_plot.mean(1,:),'Marker','.','Color', color(i,:)); hold on
% plot([1.78 2 2.22],stride_variation_speedmean_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stride_variation_speedmean_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
% for e=1:3
% stride_variation_speedmean_scatter_plot.(slopes{e})(i,:) = stride_variation_speedmean_scatter_plot.mean(e,:);
% end
% plot([0.78 1 1.22],stride_variation_speedmean_bar_plot(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([1.78 2 2.22],stride_variation_speedmean_bar_plot(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([2.78 3 3.22],stride_variation_speedmean_bar_plot(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% 
% % barplots.speedmean=bar(stride_variation_speedmean_bar_plot);
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% ylim([0.3 2])
% end 
% ylabel('Walking Speed')
% 
% 
% for i=[7 8 9 1 2 3 4 5 6]
%     if (7<=i&&i<=9)
%         err_speedvar(1,i-6)=std(variation_steps.(conds{i}).speedvar(:));
%     elseif (1<=i&&i<=3)
%         err_speedvar(2,i)=std(variation_steps.(conds{i}).speedvar(:));
%     elseif (4<=i&&i<=6)
%         err_speedvar(3,i-3)=std(variation_steps.(conds{i}).speedvar(:));
%     end 
% end
% stride_variation_speedvar_bar_plot = [mean(variation_steps.decline_050.speedvar(:)),mean(variation_steps.decline_100.speedvar(:)),mean(variation_steps.decline_150.speedvar(:));mean(variation_steps.level_050.speedvar(:)),mean(variation_steps.level_100.speedvar(:)),mean(variation_steps.level_150.speedvar(:));mean(variation_steps.incline_050.speedvar(:)),mean(variation_steps.incline_100.speedvar(:)),mean(variation_steps.incline_150.speedvar(:))];
% % barplots.speedvar=bar(stride_variation_speedvar_bar_plot);
% % set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% % hold on 
% % errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stride_variation_speedvar_bar_plot,err_speedvar,'.');
% subplot(2,4,5) 
% for i=1:s
% stride_variation_speedvar_scatter_plot.mean = [mean(variation_steps.decline_050.speedvar(i,:),2),mean(variation_steps.decline_100.speedvar(i,:),2),mean(variation_steps.decline_150.speedvar(i,:),2);mean(variation_steps.level_050.speedvar(i,:),2),mean(variation_steps.level_100.speedvar(i,:),2),mean(variation_steps.level_150.speedvar(i,:),2);mean(variation_steps.incline_050.speedvar(i,:),2),mean(variation_steps.incline_100.speedvar(i,:),2),mean(variation_steps.incline_150.speedvar(i,:),2)];
% plot([0.78 1 1.22],stride_variation_speedvar_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:)); hold on
% plot([1.78 2 2.22],stride_variation_speedvar_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stride_variation_speedvar_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
% for e=1:3
% stride_variation_speedvar_scatter_plot.(slopes{e})(i,:) = stride_variation_speedvar_scatter_plot.mean(e,:);
% end
% plot([0.78 1 1.22],stride_variation_speedvar_bar_plot(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([1.78 2 2.22],stride_variation_speedvar_bar_plot(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([2.78 3 3.22],stride_variation_speedvar_bar_plot(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% end 
% ylabel('Walking Speed Variance')

% 
% figure
% boxplot(stride_variation_speedmean_scatter_plot.(slopes{1}))
% boxplot(boxplots.speed.plot_values)
% boxplot(boxplots.speed.plot_values)
% 
% figure
% boxplots.speed.plot_values=[];
% for i=1:3
% boxplots.speed.plot_values(:,i)=[stride_variation_speedmean_scatter_plot.([slopes{i}])(:,1); stride_variation_speedmean_scatter_plot.([slopes{i}])(:,2);stride_variation_speedmean_scatter_plot.([slopes{i}])(:,3)];
% end
% % opop=boxplots.speed.plot_values;
% boxplots.speed.plot = boxplot(boxplots.speed.plot_values);
% grpstats(boxplots.speed.plot_values)
% hold on 
% for i=1:size(slopes,2)
%     for c=1:3
% %         x=repmat(i,[(size(mean(stride_variation_speedmean_scatter_plot.([slopes{i}]),1))),1]);
%         scatter(i, mean(stride_variation_speedmean_scatter_plot.([slopes{i}])(:,c)),'MarkerFaceColor',color(c*2,:)), hold on
%     end
% end

% figure 
% hh=subplot(2,2,1);
% speed_bar_plot = [mean(mean(speed_steps.decline_050(:,100:375),2)),mean(mean(speed_steps.decline_100(:,100:375),2)),mean(mean(speed_steps.decline_150(:,100:375),2));mean(mean(speed_steps.level_050(:,100:375),2)),mean(mean(speed_steps.level_100(:,100:375),2)),mean(mean(speed_steps.level_150(:,100:375),2));mean(mean(speed_steps.incline_050(:,100:375),2)),mean(mean(speed_steps.incline_100(:,100:375),2)),mean(mean(speed_steps.incline_150(:,100:375),2))];
% bar(speed_bar_plot);
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% hold on 
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],speed_bar_plot,err_speed,'.');
% for i=1:s
% speed_scatter_plot.mean = [mean(speed_steps.decline_050(i,100:375),2),mean(speed_steps.decline_100(i,100:375),2),mean(speed_steps.decline_150(i,100:375),2);mean(speed_steps.level_050(i,100:375),2),mean(speed_steps.level_100(i,100:375),2),mean(speed_steps.level_150(i,100:375),2);mean(speed_steps.incline_050(i,100:375),2),mean(speed_steps.incline_100(i,100:375),2),mean(speed_steps.incline_150(i,100:375),2)];
% % plot([0.78 1 1.22],speed_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% % plot([1.78 2 2.22],speed_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% % plot([2.78 3 3.22],speed_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
% for e=1:3
% speed_scatter_plot.(slopes{e})(i,:) = speed_scatter_plot.mean(e,:);
% end
% end 
% hh.Position=[0.3300 0.5838  0.3275 0.3412];
% ylabel('Step Speed')
% title('All Subject')

figure 
spatiotemporal_bargraphs(walking_speed,s,'walkingspeed_mean',1); 
spatiotemporal_bargraphs(walking_speed,s,'walkingspeed_var',5); 
spatiotemporal_bargraphs(frequency_steps,s,'stepfreq_mean',2); 
spatiotemporal_bargraphs(frequency_steps,s,'stepfreq_var',6); 
spatiotemporal_bargraphs(length_steps,s,'steplength_mean',3);
spatiotemporal_bargraphs(length_steps,s,'steplength_var',7);
spatiotemporal_bargraphs(width_steps,s,'stepwidth_mean',4);
spatiotemporal_bargraphs(width_steps,s,'stepwidth_var',8);

% fig3=figure; 
spatiotemporal_bargraphs_stack_controller(walking_speed,s,'walkingspeed_mean',1); 
ylim([0.8 1.8])
spatiotemporal_bargraphs_stack_controller(walking_speed,s,'walkingspeed_var',5); 
ylim([0 0.015])
spatiotemporal_bargraphs_stack_controller(frequency_steps,s,'stepfreq_mean',2); 
ylim([1.5 2.2])
% spatiotemporal_bargraphs_stack_controller(frequency_steps,s,'stepfreq_var',6); 
% ylim([0.0025 0.012])
spatiotemporal_bargraphs_stack_controller(length_steps,s,'steplength_mean',3);
ylim([0.5 0.9])
spatiotemporal_bargraphs_stack_controller(length_steps,s,'steplength_var',7);
ylim([0 0.004])
spatiotemporal_bargraphs_stack_controller(width_steps,s,'stepwidth_mean',4);
ylim([0.12 0.22])
% spatiotemporal_bargraphs_stack_controller(width_steps,s,'stepwidth_var',8);
% ylim([0.0003 0.0007])

% fig3.PaperPosition = [0.0956192152866,0.536853076490189,0.192423419043434,0.418618432199142];

% ylim([0 0.3])
% 
% gcf = fig_2;
% fig2.PaperUnits = 'inches';
% fig2.PaperPosition = [0 0 6 3];




%% NEW Figure 2 
% % hh2=figure(100);
% m=15;
% x = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5];
% figure 
% % % subplot(2,3,4)
% % stem3(x,ones([1,m])*1,mean(PowerSpectrumValues.decline_150_PowerPSD(1:m,:),2),'fill'); hold on 
% % % view([-1.419 2.2017]);
% % zlim([0 0.5])
% % xlim([0 1])
% % % xlim([0 0.1])
% % stem3(x,ones([1,m])*2,mean(PowerSpectrumValues.decline_100_PowerPSD(1:m,:),2),'fill'); hold on
% % stem3(x,ones([1,m])*3,mean(PowerSpectrumValues.decline_050_PowerPSD(1:m,:),2),'fill')
% % xlabel('Frequency')
% % % ylabel('Controller')
% % zlabel('Power')
% % title('Decline')
% 
% plot(x,mean(PowerSpectrumValues.decline_050_PowerPSD(1:m,:),2),'o-'); hold on
% plot(x,mean(PowerSpectrumValues.decline_100_PowerPSD(1:m,:),2),'o-'); hold on
% plot(x,mean(PowerSpectrumValues.decline_150_PowerPSD(1:m,:),2),'o-');  
% % view([-1.419 2.2017]);
% ylim([0 0.6])
% xlim([0 0.5])
% % xlim([0 0.1])
% 
% xlabel('Frequency')
% % ylabel('Controller')
% zlabel('Power')
% title('Decline')
% 
% 
% figure
% % subplot(2,3,5)
% plot(x,mean(PowerSpectrumValues.level_050_PowerPSD(1:m,:),2),'o-');hold on 
% plot(x,mean(PowerSpectrumValues.level_100_PowerPSD(1:m,:),2),'o-'); hold on
% plot(x,mean(PowerSpectrumValues.level_150_PowerPSD(1:m,:),2),'o-');  
% % view([-1.419 2.2017]);
% ylim([0 0.6])
% xlim([0 0.5])
% 
% % xlabel('Frequency')
% % ylabel('Controller')
% % zlabel('Power')
% title('Level')
% 
% figure
% plot(x,mean(PowerSpectrumValues.incline_050_PowerPSD(1:m,:),2),'o-');hold on 
% plot(x,mean(PowerSpectrumValues.incline_100_PowerPSD(1:m,:),2),'o-'); hold on
% plot(x,mean(PowerSpectrumValues.incline_150_PowerPSD(1:m,:),2),'o-'); 
% % view([-1.419 2.2017]);
% ylim([0 0.6])
% xlim([0 0.5])
% 
% % xlabel('Frequency')
% % ylabel('Controller')
% % zlabel('Power')
% title('Incline')
% 
% % subplot(2,3,2) 
% 
% for i=[7 8 9 1 2 3 4 5 6]
%     if (7<=i&&i<=9)
%         err_speedvar(1,i-6)=std(variation_steps.(conds{i}).speedvar(:));
%     elseif (1<=i&&i<=3)
%         err_speedvar(2,i)=std(variation_steps.(conds{i}).speedvar(:));
%     elseif (4<=i&&i<=6)
%         err_speedvar(3,i-3)=std(variation_steps.(conds{i}).speedvar(:));
%     end 
% end
% stride_variation_speedvar_bar_plot = [mean(variation_steps.decline_050.speedvar(:)),mean(variation_steps.decline_100.speedvar(:)),mean(variation_steps.decline_150.speedvar(:));mean(variation_steps.level_050.speedvar(:)),mean(variation_steps.level_100.speedvar(:)),mean(variation_steps.level_150.speedvar(:));mean(variation_steps.incline_050.speedvar(:)),mean(variation_steps.incline_100.speedvar(:)),mean(variation_steps.incline_150.speedvar(:))];
% % barplots.speedvar=bar(stride_variation_speedvar_bar_plot);
% % set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% % hold on 
% % errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stride_variation_speedvar_bar_plot,err_speedvar,'.');
% subplot(2,2,2) 
% for i=1:s
% stride_variation_speedvar_scatter_plot.mean = [mean(variation_steps.decline_050.speedvar(i,:),2),mean(variation_steps.decline_100.speedvar(i,:),2),mean(variation_steps.decline_150.speedvar(i,:),2);mean(variation_steps.level_050.speedvar(i,:),2),mean(variation_steps.level_100.speedvar(i,:),2),mean(variation_steps.level_150.speedvar(i,:),2);mean(variation_steps.incline_050.speedvar(i,:),2),mean(variation_steps.incline_100.speedvar(i,:),2),mean(variation_steps.incline_150.speedvar(i,:),2)];
% plot([0.78 1 1.22],stride_variation_speedvar_scatter_plot.mean(1,:),'Marker','.','Color',[0,0,0]+0.8); hold on
% plot([1.78 2 2.22],stride_variation_speedvar_scatter_plot.mean(2,:),'Marker','.','Color',[0,0,0]+0.8)
% plot([2.78 3 3.22],stride_variation_speedvar_scatter_plot.mean(3,:),'Marker','.','Color',[0,0,0]+0.8)
% for e=1:3
% stride_variation_speedvar_scatter_plot.(slopes{e})(i,:) = stride_variation_speedvar_scatter_plot.mean(e,:);
% end
% plot([0.78 1 1.22],stride_variation_speedvar_bar_plot(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([1.78 2 2.22],stride_variation_speedvar_bar_plot(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([2.78 3 3.22],stride_variation_speedvar_bar_plot(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% end 
% ylabel('Walking Speed Variance')
% 
% % step frequency var 
% for i=[7 8 9 1 2 3 4 5 6]
%     if (7<=i&&i<=9)
%         err_stepfreqvar(1,i-6)=std(frequency_steps.(conds{i}).stepfreq_var(:));
%     elseif (1<=i&&i<=3)
%         err_stepfreqvar(2,i)=std(frequency_steps.(conds{i}).stepfreq_var(:));
%     elseif (4<=i&&i<=6)
%         err_stepfreqvar(3,i-3)=std(frequency_steps.(conds{i}).stepfreq_var(:));
%     end 
% end
% variation_stepfreqvar_bar_plot = [mean(frequency_steps.decline_050.stepfreq_var(:)),mean(frequency_steps.decline_100.stepfreq_var(:)),mean(frequency_steps.decline_150.stepfreq_var(:));mean(frequency_steps.level_050.stepfreq_var(:)),mean(frequency_steps.level_100.stepfreq_var(:)),mean(frequency_steps.level_150.stepfreq_var(:));mean(frequency_steps.incline_050.stepfreq_var(:)),mean(frequency_steps.incline_100.stepfreq_var(:)),mean(frequency_steps.incline_150.stepfreq_var(:))];
% % barplots.stepfreq_var=bar(variation_stepfreqvar_bar_plot);
% % set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% % hold on 
% % errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],variation_stepfreqvar_bar_plot,err_stepfreqvar,'.');
% subplot(2,2,1) 
% for i=1:s
% variation_stepfreqvar_scatter_plot.mean = [mean(frequency_steps.decline_050.stepfreq_var(i,:),2),mean(frequency_steps.decline_100.stepfreq_var(i,:),2),mean(frequency_steps.decline_150.stepfreq_var(i,:),2);mean(frequency_steps.level_050.stepfreq_var(i,:),2),mean(frequency_steps.level_100.stepfreq_var(i,:),2),mean(frequency_steps.level_150.stepfreq_var(i,:),2);mean(frequency_steps.incline_050.stepfreq_var(i,:),2),mean(frequency_steps.incline_100.stepfreq_var(i,:),2),mean(frequency_steps.incline_150.stepfreq_var(i,:),2)];
% plot([0.78 1 1.22],variation_stepfreqvar_scatter_plot.mean(1,:),'Marker','.','Color',[0,0,0]+0.8); hold on
% plot([1.78 2 2.22],variation_stepfreqvar_scatter_plot.mean(2,:),'Marker','.','Color',[0,0,0]+0.8)
% plot([2.78 3 3.22],variation_stepfreqvar_scatter_plot.mean(3,:),'Marker','.','Color',[0,0,0]+0.8)
% for e=1:3
% variation_stepfreqvar_scatter_plot.(slopes{e})(i,:) = variation_stepfreqvar_scatter_plot.mean(e,:);
% end
% plot([0.78 1 1.22],variation_stepfreqvar_bar_plot(1,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([1.78 2 2.22],variation_stepfreqvar_bar_plot(2,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% plot([2.78 3 3.22],variation_stepfreqvar_bar_plot(3,:),'linewidth',2,'Marker','.','MarkerSize',12,'Color',[0,0,0])
% set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'});
% end 
% ylabel('Step Frequency Variance')

%% NEW Figure 2.5 - box
figure
boxplots.steplengthvar_speedtrend.plot_values=[];
for i=1
boxplots.steplengthvar_speedtrend.plot_values=[stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,1) stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,2) stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,3)];
end
% opop=boxplots.speed.plot_values;
boxplots.steplengthvar_speedtrend.plot = boxplot(boxplots.steplengthvar_speedtrend.plot_values,'Whisker',3);
hold on 
for i=1:size(slopes,2)
    for c=1:3
%         x=repmat(i,[(size(mean(stride_variation_speedmean_scatter_plot.([slopes{i}]),1))),1]);
        scatter(i,stridelength_variation_speedtrend_bar_plot(i,c),'MarkerFaceColor',color(c*2,:)), hold on
    end
end 

figure
boxplots.steplengthvar_speedtrend.plot_values=[];
for i=1:3
boxplots.steplengthvar_speedtrend.plot_values(:,i)=[stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,1); stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,2);stridelength_variation_speedtrend_scatter_plot.([slopes{i}])(:,3)];
end
% opop=boxplots.speed.plot_values;
boxplots.steplengthvar_speedtrend.plot = boxplot(boxplots.steplengthvar_speedtrend.plot_values);
hold on 
for i=1:size(slopes,2)
    for c=1:3
%         x=repmat(i,[(size(mean(stride_variation_speedmean_scatter_plot.([slopes{i}]),1))),1]);
        scatter(i,stridelength_variation_speedtrend_bar_plot(i,c),'MarkerFaceColor',color(c*2,:)), hold on
    end
end 

%% NEW Figure 2.5 
stridelength_stacked=[];
for i=[1 4 7 2 5 8 3 6 9]
e=size(stridelength_stacked,1)+1;
stridelength_stacked(e,:)=[stridelength_variation_detrended_bar_plot(i),stridelength_variation_speedtrend_bar_plot(i)];
stridelength_stacked(e+1,:)=[stridelength_variation_total_bar_plot(i),0]; 
end 
figure;
hh3=subplot(2,2,1);
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
stridelength_variation_total_scatter_plot.mean = [mean(frequency_steps.decline_050.totalvar(i,:),2),mean(variation_steps.decline_100.totalvar(i,:),2),mean(variation_steps.decline_150.totalvar(i,:),2);mean(variation_steps.level_050.totalvar(i,:),2),mean(variation_steps.level_100.totalvar(i,:),2),mean(variation_steps.level_150.totalvar(i,:),2);mean(variation_steps.incline_050.totalvar(i,:),2),mean(variation_steps.incline_100.totalvar(i,:),2),mean(variation_steps.incline_150.totalvar(i,:),2)];
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
hh3.Position=[0.1386 0.5838 0.7661 0.3390];
ylim([0 0.004])
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Length Variability')


stridewidth_stacked=[];
for i=[1 4 7 2 5 8 3 6 9]
e=size(stridewidth_stacked,1)+1;
stridewidth_stacked(e,:)=[stridewidth_variation_detrended_bar_plot(i),stridewidth_variation_speedtrend_bar_plot(i)];
stridewidth_stacked(e+1,:)=[stridewidth_variation_total_bar_plot(i),0]; 
end 
% figure;
hh4=subplot(2,2,3);
h=bar([0.80 0.90 1.15 1.25 1.50 1.60],stridewidth_stacked(1:6,:),'stacked'); hold on
colors_stacks = jet(size(h,2));
colors_stacks = repelem(colors_stacks,size(h,1),1); 
colors_stacks = mat2cell(colors_stacks,ones(size(colors_stacks,1),1),3);
set(h,{'FaceColor'},colors_stacks)
h=bar([2.60 2.70 2.95 3.05 3.30 3.40],stridewidth_stacked(7:12,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
hold on
h=bar([4.40 4.50 4.75 4.85 5.10 5.20],stridewidth_stacked(13:18,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
set(gca, 'XTick', [1.2,3,4.8],'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
for i=1:s
stridewidth_variation_total_scatter_plot.mean = [mean(variation_steps_width.decline_050.totalvar(i,:),2),mean(variation_steps_width.decline_100.totalvar(i,:),2),mean(variation_steps_width.decline_150.totalvar(i,:),2);mean(variation_steps_width.level_050.totalvar(i,:),2),mean(variation_steps_width.level_100.totalvar(i,:),2),mean(variation_steps_width.level_150.totalvar(i,:),2);mean(variation_steps_width.incline_050.totalvar(i,:),2),mean(variation_steps_width.incline_100.totalvar(i,:),2),mean(variation_steps_width.incline_150.totalvar(i,:),2)];
% plot([0.90 1.25 1.60],stridewidth_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([2.70 3.05 3.40],stridewidth_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([4.50 4.85 5.20],stridewidth_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_total_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_total_scatter_plot.mean(e,:);
end
end 
hold on 
errorbar([0.90 1.25 1.60; 2.70 3.05 3.40; 4.50 4.85 5.20],stridewidth_variation_total_bar_plot,error_width_totalvar,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
hh4.Position=[0.1386 0.1100 0.7661 0.3390];
ylim([0 0.0008])
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Width Variability')

%% NEW Figure - Controller on x-axis
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_detrended(1,i-6)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (1<=i&&i<=3)
        error_length_detrended(2,i)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (4<=i&&i<=6)
        error_length_detrended(3,i-3)=std(variation_steps.(conds{i}).slminusfit(:));
    end 
end
stridelength_stacked=[];
for i=1:9
e=size(stridelength_stacked,1)+1;
stridelength_stacked(e,:)=[stridelength_variation_detrended_bar_plot(i),stridelength_variation_speedtrend_bar_plot(i)];
stridelength_stacked(e+1,:)=[stridelength_variation_total_bar_plot(i),0]; 
end 
figure;
hh3=subplot(2,2,1);
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
set(gca, 'XTick', [1.2,3,4.8],'XTickLabel',{'Low' 'Medium' 'High'});
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
errorbar([0.90 1.25 1.60; 2.70 3.05 3.40; 4.50 4.85 5.20],stridelength_variation_total_bar_plot',error_length_totalvar','.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridelength_variation_detrended_bar_plot',error_length_detrended','.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridelength_variation_detrended_bar_plot' + stridelength_variation_speedtrend_bar_plot'),error_length_speedtrend','.'); 
hh3.Position=[0.1386 0.5838 0.7661 0.3390];
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Length Variability')


stridewidth_stacked=[];
for i=1:9
e=size(stridewidth_stacked,1)+1;
stridewidth_stacked(e,:)=[stridewidth_variation_detrended_bar_plot(i),stridewidth_variation_speedtrend_bar_plot(i)];
stridewidth_stacked(e+1,:)=[stridewidth_variation_total_bar_plot(i),0]; 
end 
% figure;
hh4=subplot(2,2,3);
h=bar([0.80 0.90 1.15 1.25 1.50 1.60],stridewidth_stacked(1:6,:),'stacked'); hold on
colors_stacks = jet(size(h,2));
colors_stacks = repelem(colors_stacks,size(h,1),1); 
colors_stacks = mat2cell(colors_stacks,ones(size(colors_stacks,1),1),3);
set(h,{'FaceColor'},colors_stacks)
h=bar([2.60 2.70 2.95 3.05 3.30 3.40],stridewidth_stacked(7:12,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
hold on
h=bar([4.40 4.50 4.75 4.85 5.10 5.20],stridewidth_stacked(13:18,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
set(gca, 'XTick', [1.2,3,4.8],'XTickLabel',{'Low' 'Medium' 'High'});
hold on 
for i=1:s
stridewidth_variation_total_scatter_plot.mean = [mean(variation_steps_width.decline_050.totalvar(i,:),2),mean(variation_steps_width.decline_100.totalvar(i,:),2),mean(variation_steps_width.decline_150.totalvar(i,:),2);mean(variation_steps_width.level_050.totalvar(i,:),2),mean(variation_steps_width.level_100.totalvar(i,:),2),mean(variation_steps_width.level_150.totalvar(i,:),2);mean(variation_steps_width.incline_050.totalvar(i,:),2),mean(variation_steps_width.incline_100.totalvar(i,:),2),mean(variation_steps_width.incline_150.totalvar(i,:),2)];
% plot([0.90 1.25 1.60],stridewidth_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([2.70 3.05 3.40],stridewidth_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([4.50 4.85 5.20],stridewidth_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_total_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_total_scatter_plot.mean(e,:);
end
end 
hold on 
errorbar([0.90 1.25 1.60; 2.70 3.05 3.40; 4.50 4.85 5.20],stridewidth_variation_total_bar_plot',error_width_totalvar','.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot',error_width_detrended','.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot' + stridewidth_variation_speedtrend_bar_plot'),error_width_speedtrend','.');
hh4.Position=[0.1386 0.1100 0.7661 0.3390];
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Width Variability')


%% NEW Figure 2.6

for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_detrended(1,i-6)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (1<=i&&i<=3)
        error_length_detrended(2,i)=std(variation_steps.(conds{i}).slminusfit(:));
    elseif (4<=i&&i<=6)
        error_length_detrended(3,i-3)=std(variation_steps.(conds{i}).slminusfit(:));
    end 
end
stridelength_stacked=[];
for i=[1 4 7 2 5 8 3 6 9]
e=size(stridelength_stacked,1)+1;
stridelength_stacked(e,:)=[stridelength_variation_detrended_bar_plot(i),stridelength_variation_speedtrend_bar_plot(i)];
stridelength_stacked(e+1,:)=[stridelength_variation_total_bar_plot(i),0]; 
end 
figure;
hh3=subplot(2,2,1);
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
hh3.Position=[0.1386 0.5838 0.7661 0.3390];
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Length Variability')

subplot(2,2,3)
bar(stridelength_variation_detrended_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'}); hold on,
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_detrended_bar_plot,error_length_speedtrend,'.');
ylabel('Detrended Step Length')
ylim([0 0.0015])

subplot(2,2,4)
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_length_speedtrend(1,i-6)=std(variation_steps.(conds{i}).speedtrend(:));
    elseif (1<=i&&i<=3)
        error_length_speedtrend(2,i)=std(variation_steps.(conds{i}).speedtrend(:));
    elseif (4<=i&&i<=6)
        error_length_speedtrend(3,i-3)=std(variation_steps.(conds{i}).speedtrend(:));
    end 
end
bar(stridelength_variation_speedtrend_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'}); hold on,
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_speedtrend_bar_plot,error_length_speedtrend,'.');
ylabel('Speedtrend Step Length')
ylim([0 0.0015])

 %% NEW Figure 2.7
stridewidth_stacked=[];
for i=[1 4 7 2 5 8 3 6 9]
e=size(stridewidth_stacked,1)+1;
stridewidth_stacked(e,:)=[stridewidth_variation_detrended_bar_plot(i),stridewidth_variation_speedtrend_bar_plot(i)];
stridewidth_stacked(e+1,:)=[stridewidth_variation_total_bar_plot(i),0]; 
end 
figure;
hh4=subplot(2,2,1);
h=bar([0.80 0.90 1.15 1.25 1.50 1.60],stridewidth_stacked(1:6,:),'stacked'); hold on
colors_stacks = jet(size(h,2));
colors_stacks = repelem(colors_stacks,size(h,1),1); 
colors_stacks = mat2cell(colors_stacks,ones(size(colors_stacks,1),1),3);
set(h,{'FaceColor'},colors_stacks)
h=bar([2.60 2.70 2.95 3.05 3.30 3.40],stridewidth_stacked(7:12,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
hold on
h=bar([4.40 4.50 4.75 4.85 5.10 5.20],stridewidth_stacked(13:18,:),'stacked');
set(h,{'FaceColor'},colors_stacks)
set(gca, 'XTick', [1.2,3,4.8],'XTickLabel',{'Decline' 'Level' 'Incline'});
hold on 
for i=1:s
stridewidth_variation_total_scatter_plot.mean = [mean(variation_steps_width.decline_050.totalvar(i,:),2),mean(variation_steps_width.decline_100.totalvar(i,:),2),mean(variation_steps_width.decline_150.totalvar(i,:),2);mean(variation_steps_width.level_050.totalvar(i,:),2),mean(variation_steps_width.level_100.totalvar(i,:),2),mean(variation_steps_width.level_150.totalvar(i,:),2);mean(variation_steps_width.incline_050.totalvar(i,:),2),mean(variation_steps_width.incline_100.totalvar(i,:),2),mean(variation_steps_width.incline_150.totalvar(i,:),2)];
% plot([0.90 1.25 1.60],stridewidth_variation_total_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([2.70 3.05 3.40],stridewidth_variation_total_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([4.50 4.85 5.20],stridewidth_variation_total_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_total_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_total_scatter_plot.mean(e,:);
end
end 
hold on 
errorbar([0.90 1.25 1.60; 2.70 3.05 3.40; 4.50 4.85 5.20],stridewidth_variation_total_bar_plot,error_width_totalvar,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],stridewidth_variation_detrended_bar_plot,error_width_detrended,'.');
errorbar([0.80 1.15 1.50; 2.60 2.95 3.30; 4.40 4.75 5.10],(stridewidth_variation_detrended_bar_plot + stridewidth_variation_speedtrend_bar_plot),error_width_speedtrend,'.');
hh4.Position=[0.1386 0.5838 0.7661 0.3390];
ylabel('Variance (m^2)')
legend('Detrended','Speedtrend')
title('Step Width Variability')

subplot(2,2,3)
bar(stridewidth_variation_detrended_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'}); hold on,
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_detrended_bar_plot,error_length_speedtrend,'.');
ylabel('Detrended Step Width')
ylim([0 0.0006])

subplot(2,2,4)
bar(stridewidth_variation_speedtrend_bar_plot)
set(gca, 'XTick', 1:3,'XTickLabel',{'Decline' 'Level' 'Incline'}); hold on,
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_speedtrend_bar_plot,error_length_speedtrend,'.');
ylabel('Speedtrend Step Width')
ylim([0 0.0006])

%% NEW Figure 3 
figure
subplot(1,2,1) 
stridelength_variation_detrended_vs_slopes = [mean(variation_steps.decline_050.slminusfit(:)),mean(variation_steps.level_050.slminusfit(:)),mean(variation_steps.incline_050.slminusfit(:));mean(variation_steps.decline_100.slminusfit(:)),mean(variation_steps.level_100.slminusfit(:)),mean(variation_steps.incline_100.slminusfit(:));mean(variation_steps.decline_150.slminusfit(:)),mean(variation_steps.level_150.slminusfit(:)),mean(variation_steps.incline_150.slminusfit(:))];
bar(stridelength_variation_detrended_vs_slopes)
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_detrended_vs_slopes,error_length_detrended_fig,'.');
for i=1:s
stridelength_variation_detrended_vs_slope_scatter_plot.mean = [mean(variation_steps.decline_050.slminusfit(i,:),2),mean(variation_steps.level_050.slminusfit(i,:),2),mean(variation_steps.incline_050.slminusfit(i,:),2);mean(variation_steps.decline_100.slminusfit(i,:),2),mean(variation_steps.level_100.slminusfit(i,:),2),mean(variation_steps.incline_100.slminusfit(i,:),2);mean(variation_steps.decline_150.slminusfit(i,:),2),mean(variation_steps.level_150.slminusfit(i,:),2),mean(variation_steps.incline_150.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_detrended_vs_slope_scatter_plot.(slopes{e})(i,:) = stridelength_variation_detrended_vs_slope_scatter_plot.mean(e,:);
end
end 
ylabel('Step Length Variability -- Detrended vs Slopes')
ylim([0 0.002])
legend('Decline','Level','Incline')

subplot(1,2,2) 
stridelength_variation_speedtrend_vs_slopes = [mean(variation_steps.decline_050.speedtrend(:)),mean(variation_steps.level_050.speedtrend(:)),mean(variation_steps.incline_050.speedtrend(:));mean(variation_steps.decline_100.speedtrend(:)),mean(variation_steps.level_100.speedtrend(:)),mean(variation_steps.incline_100.speedtrend(:));mean(variation_steps.decline_150.speedtrend(:)),mean(variation_steps.level_150.speedtrend(:)),mean(variation_steps.incline_150.speedtrend(:))];
bar(stridelength_variation_speedtrend_vs_slopes)
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
% errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridelength_variation_speedtrend_vs_slopes,error_length_detrended_fig,'.');
for i=1:s
stridelength_variation_speedtrend_vs_slope_scatter_plot.mean = [mean(variation_steps.decline_050.speedtrend(i,:),2),mean(variation_steps.level_050.speedtrend(i,:),2),mean(variation_steps.incline_050.speedtrend(i,:),2);mean(variation_steps.decline_100.speedtrend(i,:),2),mean(variation_steps.level_100.speedtrend(i,:),2),mean(variation_steps.incline_100.speedtrend(i,:),2);mean(variation_steps.decline_150.speedtrend(i,:),2),mean(variation_steps.level_150.speedtrend(i,:),2),mean(variation_steps.incline_150.speedtrend(i,:),2)];
% plot([0.78 1 1.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridelength_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridelength_variation_speedtrend_vs_slope_scatter_plot.(slopes{e})(i,:) = stridelength_variation_speedtrend_vs_slope_scatter_plot.mean(e,:);
end
end 
ylabel('Step Length Variability -- Speedtrend vs Slopes')
ylim([0 0.002])
legend('Decline','Level','Incline')

%% NEW Figure 4 
figure
subplot(1,2,1)
jj=[7 1 4 8 2 5 9 3 6];
for i=[7 1 4 8 2 5 9 3 6]
    m=find(jj==i);
    if m<4
        error_width_detrended_fig(1,m)=error_width_detrended(i);
    elseif m<7
        error_width_detrended_fig(2,m-3)=error_width_detrended(i);
    else
        error_width_detrended_fig(3,m-6)=error_width_detrended(i);
    end
end

stridewidth_variation_detrended_vs_slopes = [mean(variation_steps_width.decline_050.slminusfit(:)),mean(variation_steps_width.level_050.slminusfit(:)),mean(variation_steps_width.incline_050.slminusfit(:));mean(variation_steps_width.decline_100.slminusfit(:)),mean(variation_steps_width.level_100.slminusfit(:)),mean(variation_steps_width.incline_100.slminusfit(:));mean(variation_steps_width.decline_150.slminusfit(:)),mean(variation_steps_width.level_150.slminusfit(:)),mean(variation_steps_width.incline_150.slminusfit(:))];
bar(stridewidth_variation_detrended_vs_slopes)
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridewidth_variation_detrended_vs_slopes,error_width_detrended_fig,'.');
for i=1:s
stridewidth_variation_detrended_vs_slope_scatter_plot.mean = [mean(variation_steps_width.decline_050.slminusfit(i,:),2),mean(variation_steps_width.level_050.slminusfit(i,:),2),mean(variation_steps_width.incline_050.slminusfit(i,:),2);mean(variation_steps_width.decline_100.slminusfit(i,:),2),mean(variation_steps_width.level_100.slminusfit(i,:),2),mean(variation_steps_width.incline_100.slminusfit(i,:),2);mean(variation_steps_width.decline_150.slminusfit(i,:),2),mean(variation_steps_width.level_150.slminusfit(i,:),2),mean(variation_steps_width.incline_150.slminusfit(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_detrended_vs_slope_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_detrended_vs_slope_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Detrended vs Slopes')
ylim([0 0.0008])
legend('Decline','Level','Incline')

subplot(1,2,2)
for i=[7 8 9 1 2 3 4 5 6]
    if (7<=i&&i<=9)
        error_width_speedtrend(1,i-6)=std(variation_steps_width.(conds{i}).speedtrend(:));
    elseif (1<=i&&i<=3)
        error_width_speedtrend(2,i)=std(variation_steps_width.(conds{i}).speedtrend(:));
    elseif (4<=i&&i<=6)
        error_width_speedtrend(3,i-3)=std(variation_steps_width.(conds{i}).speedtrend(:));
    end 
end
stridewidth_variation_speedtrend_vs_slopes = [mean(variation_steps_width.decline_050.speedtrend(:)),mean(variation_steps_width.level_050.speedtrend(:)),mean(variation_steps_width.incline_050.speedtrend(:));mean(variation_steps_width.decline_100.speedtrend(:)),mean(variation_steps_width.level_100.speedtrend(:)),mean(variation_steps_width.incline_100.speedtrend(:));mean(variation_steps_width.decline_150.speedtrend(:)),mean(variation_steps_width.level_150.speedtrend(:)),mean(variation_steps_width.incline_150.speedtrend(:))];
bar(stridewidth_variation_speedtrend_vs_slopes)
set(gca, 'XTick', 1:3,'XTickLabel',{'0.50' '1.00' '1.50'});
hold on 
errorbar([0.78 1 1.22; 1.78 2 2.22; 2.78 3 3.22],stridewidth_variation_speedtrend_vs_slopes,error_width_speedtrend,'.');
for i=1:s
stridewidth_variation_speedtrend_vs_slope_scatter_plot.mean = [mean(variation_steps_width.decline_050.speedtrend(i,:),2),mean(variation_steps_width.level_050.speedtrend(i,:),2),mean(variation_steps_width.incline_050.speedtrend(i,:),2);mean(variation_steps_width.decline_100.speedtrend(i,:),2),mean(variation_steps_width.level_100.speedtrend(i,:),2),mean(variation_steps_width.incline_100.speedtrend(i,:),2);mean(variation_steps_width.decline_150.speedtrend(i,:),2),mean(variation_steps_width.level_150.speedtrend(i,:),2),mean(variation_steps_width.incline_150.speedtrend(i,:),2)];
% plot([0.78 1 1.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(1,:),'Marker','.','Color',color(i,:))
% plot([1.78 2 2.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(2,:),'Marker','.','Color',color(i,:))
% plot([2.78 3 3.22],stridewidth_variation_detrended_vs_slope_scatter_plot.mean(3,:),'Marker','.','Color',color(i,:))
for e=1:3
stridewidth_variation_speedtrend_vs_slope_scatter_plot.(slopes{e})(i,:) = stridewidth_variation_speedtrend_vs_slope_scatter_plot.mean(e,:);
end
end 
ylabel('Step Width Variability -- Speedtrend vs Slopes')
ylim([0 0.0008])
legend('Decline','Level','Incline')


% speed_bar_plot = [mean(mean(speed_steps.decline_050(:,100:375),2)),mean(mean(speed_steps.decline_100(:,100:375),2)),mean(mean(speed_steps.decline_150(:,100:375),2));mean(mean(speed_steps.level_050(:,100:375),2)),mean(mean(speed_steps.level_100(:,100:375),2)),mean(mean(speed_steps.level_150(:,100:375),2));mean(mean(speed_steps.incline_050(:,100:375),2)),mean(mean(speed_steps.incline_100(:,100:375),2)),mean(mean(speed_steps.incline_150(:,100:375),2))];
% variation_steps.(conds{c}).totalstd(s,:) 

%% single subject plots 

%% Stats 
% stats_names = {'speed' 'speedvar' 'stridelength_detrended' 'stridelength_totalvar' 'stridelength_speedtrend' 'stridelength_slopes_detrended' 'stridewidth_slopes_detrended' 'stridewidth_speedtrend' 'stridewidth_detrended' 'stridewidth_totalvar'};
% for i=1:length(stats_names)
%         stats.(stats_names{i}).allslopes=table('Size',[s*3 3],'VariableTypes',repmat("double",1,3),'VariableNames',["Condition_050" "Condition_100" "Condition_150"]);
% end
% 
% new=[1,11,21];
% for m=1:3
%         for e=1:3       
%             stats.speed.allsensitivity_matrix(new(m):new(m)+9,e)=stride_variation_speedmean_scatter_plot.(slopes{m})(:,e); 
% %             stats.speedvar.allslopes{i,e}=stride_variation_speedvar_scatter_plot.(slopes{m});
% %             stats.stridelength_detrended.allslopes{i,e}=stridelength_variation_detrended_scatter_plot.(slopes{m}); 
% %             stats.stridelength_totalvar.allslopes{i,e}=stridelength_variation_total_scatter_plot.(slopes{m});
% %             stats.stridelength_speedtrend.allslopes{i,e}=stridelength_variation_speedtrend_scatter_plot.(slopes{m});
% % %             stats.stridelength_slopes_detrended.allslopes{i,e}=stridelength_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% % %             stats.stridewidth_slopes_detrended.allslopes{i,e}=stridewidth_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% %             stats.stridewidth_speedtrend.allslopes{i,e}=stridewidth_variation_speedtrend_scatter_plot.(slopes{m});
% %             stats.stridewidth_totalvar.allslopes{i,e}=stridewidth_variation_total_scatter_plot.(slopes{m});
% %             stats.stridewidth_detrended.allslopes{i,e}=stridewidth_variation_detrended_scatter_plot.(slopes{m});
%         end
% end 
%   
% for m=1:3     
%     stats.speed.allslopes{:,m}= stats.speed.allsensitivity_matrix(:,m);
% %             stats.speedvar.allslopes{i,e}=stride_variation_speedvar_scatter_plot.(slopes{m});
% %             stats.stridelength_detrended.allslopes{i,e}=stridelength_variation_detrended_scatter_plot.(slopes{m}); 
% %             stats.stridelength_totalvar.allslopes{i,e}=stridelength_variation_total_scatter_plot.(slopes{m});
% %             stats.stridelength_speedtrend.allslopes{i,e}=stridelength_variation_speedtrend_scatter_plot.(slopes{m});
% % %             stats.stridelength_slopes_detrended.allslopes{i,e}=stridelength_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% % %             stats.stridewidth_slopes_detrended.allslopes{i,e}=stridewidth_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% %             stats.stridewidth_speedtrend.allslopes{i,e}=stridewidth_variation_speedtrend_scatter_plot.(slopes{m});
% %             stats.stridewidth_totalvar.allslopes{i,e}=stridewidth_variation_total_scatter_plot.(slopes{m});
% %             stats.stridewidth_detrended.allslopes{i,e}=stridewidth_variation_detrended_scatter_plot.(slopes{m});
% 
% end 
%      
% for m=1:3     
%     stats.speed.allslopes.Slopes(new(m):new(m)+9,1)=categorical(cellstr(repmat(slopes(m),10,1)));
% %             stats.speedvar.allslopes{i,e}=stride_variation_speedvar_scatter_plot.(slopes{m});
% %             stats.stridelength_detrended.allslopes{i,e}=stridelength_variation_detrended_scatter_plot.(slopes{m}); 
% %             stats.stridelength_totalvar.allslopes{i,e}=stridelength_variation_total_scatter_plot.(slopes{m});
% %             stats.stridelength_speedtrend.allslopes{i,e}=stridelength_variation_speedtrend_scatter_plot.(slopes{m});
% % %             stats.stridelength_slopes_detrended.allslopes{i,e}=stridelength_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% % %             stats.stridewidth_slopes_detrended.allslopes{i,e}=stridewidth_variation_detrended_vs_slope_scatter_plot.(slopes{m});
% %             stats.stridewidth_speedtrend.allslopes{i,e}=stridewidth_variation_speedtrend_scatter_plot.(slopes{m});
% %             stats.stridewidth_totalvar.allslopes{i,e}=stridewidth_variation_total_scatter_plot.(slopes{m});
% %             stats.stridewidth_detrended.allslopes{i,e}=stridewidth_variation_detrended_scatter_plot.(slopes{m});
% 
% end 
%   
% for i=1:length(stats_names)
%         stats.(stats_names{i}).rm=fitrm([stats.(stats_names{i}).allslopes],'measure~Slopes');
%         stats.(stats_names{i}).rmResult=ranova(stats.(stats_names{i}).rm);
%         stats.(stats_names{i}).pairtest=multcompare(stats.(stats_names{i}).rm,'Time','ComparisonType','tukey-kramer');
% end
% 
% % compare slopes 
% stats.speed_allslopes.slopes_means = table('Size',[s 3],'VariableTypes',repmat("double",1,3),'VariableNames',["Decline" "Level" "Incline"]);
% for e=1:3 
% stats.speed_allslopes.slopes_means{:,e} = mean(stride_variation_speedmean_scatter_plot.(slopes{e}),2);
% end 
% for e=1:3
%         stats.speed_allslopes.rm.slopes_means=fitrm(stats.speed_allslopes.slopes_means,'Decline-Incline~1');
%         stats.speed_allslopes.rmResult.slopes_means=ranova(stats.speed_allslopes.rm.slopes_means);
%         stats.speed_allslopes.pairtest.slopes_means=multcompare(stats.speed_allslopes.rm.slopes_means,'Time','ComparisonType','tukey-kramer');
% end
%% plot to explain fit 

[curvefit,gof,output] = fit(cdate,pop,'poly3','normalize','on')

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
ylim([0.5 1.2])
xlim([0 325])
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
ylim([0.4 0.7])
xlim([0.7 1.3])

% figure
subplot(4,4,8+i)
plot(fitplot_actual_steplength.(conds{i}).(subjs{m}),'.','Color','k'),hold on
plot(fitplot_fitted_steplength.(conds{i}).(subjs{m}),'O','Color','g','LineWidth',2), hold on  

if i==3
legend('Actual Step length','Speed-trend')
end
xlim([0 325])
ylabel('Step Length (m)')
xlabel('Steps')
ylim([0.4 0.7])


% figure
subplot(4,4,12+i)
plot(fitplot_stepfitted_minus_actualstep.(conds{i}).(subjs{m}),'.','Color','g','MarkerSize',12), hold on
xlim([0 325])
ylabel('Detrended')
xlabel('Steps')
ylim([-0.075 0.075])

end
end

for i=1:3
% figure(30) 
% figure1 = figure('WindowState','maximized');
box(subplot(4,4,8+i),'on');
txt_on_value2=num2str(variation_steps.(conds{i}).speedtrend(m,:));
merge_text2={'Speed-trend variance',txt_on_value2};
annotation(figure(30),'textbox',...
    [0.004125+(i*0.2) 0.330456905503634 0.0495208333333333 0.0218068535825545],...
    'String',merge_text2,...
    'FitBoxToText','on');
end

for i=1:3
    

% figure(30) 
% figure1 = figure('WindowState','maximized');
box(subplot(4,4,12+i),'on');
txt_on_value=num2str(variation_steps.(conds{i}).slminusfit(m,:));
merge_text={'Detrended variance',txt_on_value};
annotation(figure(30),'textbox',[0.004125+(i*0.2) 0.12456905503634 0.0495208333333333 0.0218068535825545],...
    'String',merge_text,...
    'FitBoxToText','on');
end


% [counts,bins] = hist(testData); %# get counts and bin locations
% barh(bins,counts)
% saveas(figure(15),'explain_detrend.fig')

%% github changes 
% Stats 
stats_names = {'speed' 'speedvar' 'stridelength_detrended' 'stridelength_totalvar' 'stridelength_speedtrend' 'stridelength_slopes_detrended' 'stridewidth_slopes_detrended' 'stridewidth_speedtrend' 'stridewidth_detrended' 'stridewidth_totalvar'};
for i=1:length(stats_names)
    for e=1:3
        stats.(stats_names{i}).(slopes{e})=table('Size',[s 3],'VariableTypes',repmat("double",1,3),'VariableNames',["Condition_050" "Condition_100" "Condition_150"]);
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
        stats.(stats_names{i}).pairtest.(slopes{e})=multcompare(stats.(stats_names{i}).rm.(slopes{e}),'Time','ComparisonType','tukey-kramer');
    end
end

% compare slopes 
stats.speed_allslopes.slopes_means = table('Size',[s 3],'VariableTypes',repmat("double",1,3),'VariableNames',["Decline" "Level" "Incline"]);
for e=1:3 
stats.speed_allslopes.slopes_means{:,e} = mean(stride_variation_speedmean_scatter_plot.(slopes{e}),2);
end 
for e=1:3
        stats.speed_allslopes.rm.slopes_means=fitrm(stats.speed_allslopes.slopes_means,'Decline-Incline~1');
        stats.speed_allslopes.rmResult.slopes_means=ranova(stats.speed_allslopes.rm.slopes_means);
        stats.speed_allslopes.pairtest.slopes_means=multcompare(stats.speed_allslopes.rm.slopes_means,'Time','ComparisonType','tukey-kramer');
end
%% for minitab
% for i=1:s
%     for m=1:c
%         minitab.stepfreq_mean = frequency_steps.(conds{m}).stepfreq_var(i);
%     end
% end

        