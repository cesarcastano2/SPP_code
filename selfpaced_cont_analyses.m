clearvars; clc
close all

g = 1;
P = -1;
D = 1;
% ae = [0.2 0.999];
% be = [-0.001 -0.0001];

ae = [1 0.999];
be = [-0.0001 -0.0001];

normmag = 1; 

fs = 240; %hz, resamp tm to df rate


subjs = {'s20'};
conds = {'no_perturb' 'no_mf' 'diff_f' 'diff_m' 'diff_fm'};
projfolder='F:\SPP\Helen'; %for new data

% condlabel = ["walk", "stand, start back", "stand, start middle", "stand, start front"];

% gain = [1 1.5 0.5];

% figure(1)
% set(1,'color','white');
% 
% figure(2)
% set(2,'color','white','name',"eqn terms");

for s = 1:length(subjs)
    for c = 1:length(conds) 
        dflow_file = [projfolder '/' subjs{s} '_' conds{c} '0001.txt'];
        dflow_treadmill_file = [projfolder '/' subjs{s} '_' conds{c} '_treadmill0001.txt'];
        
        tm.(conds{c}) = importTreadmillFile_SPP(dflow_file);
        df.(conds{c}) = importDflowFile(dflow_file);
    end
end


for c = 1:13
    
    filename = "subj" + "_" + "condition" + num2str(c) + "_treadmill0010001.txt";
    tm = importTreadmillFile(dflow_file);
    
    filename = "condition" + num2str(c) + "0001.txt";
    df = importDflowFile(dflow_file);

    DATA{c}.belt.time = interp1(tm{:,1},tm{:,1},df{:,1},'spline');
    DATA{c}.belt.speed.raw = interp1(tm{:,1},tm{:,2},df{:,1},'spline');
    
    istart = find(abs(DATA{c}.belt.speed.raw) > 0.01,1,'first');
    
    DATA{c}.belt.time(1:istart-1) = [];
    DATA{c}.belt.speed.raw(1:istart-1) = [];
    
    DATA{c}.belt.time = DATA{c}.belt.time-DATA{c}.belt.time(1);

    for m = ["RASI" "LASI" "RPSI" "LPSI" "RHEE"]
        M.(m) = [df{istart:end,m + "PosX"}, df{istart:end,m + "PosY"}, df{istart:end,m + "PosZ"}];
    end

    DATA{c}.com.raw = mean(cat(3,M.RASI, M.LASI, M.RPSI, M.LPSI),3);
    
    fc = 6; % low pass freq cutoff
    [b,a] = butter(2,fc/(fs/2));
    DATA{c}.com.lpf6 = filtfilt(b,a,DATA{c}.com.raw); 
    DATA{c}.rheel = filtfilt(b,a,M.RHEE);
    
    DATA{c}.time = df{istart:end,"TimeStamp"}-df{istart,"TimeStamp"};
    
    fc = 0.5; % low pass freq cutoff
    [b,a] = butter(2,fc/(fs/2));
    DATA{c}.com.lpf05 = filtfilt(b,a,DATA{c}.com.lpf6);
    DATA{c}.com.dot = diff(DATA{c}.com.lpf05)*fs;
    DATA{c}.com.ddot = diff(DATA{c}.com.dot)*fs;

    DATA{c}.belt.speed.lpf05 = filtfilt(b,a,DATA{c}.belt.speed.raw);
    DATA{c}.belt.speed.dot = filtfilt(b,a,diff(DATA{c}.belt.speed.lpf05)*fs);
    DATA{c}.belt.speed.ddot = filtfilt(b,a,diff(DATA{c}.belt.speed.dot)*fs);
    
    DATA{c}.correction_PxD = g*(P*DATA{c}.com.lpf05(1:end-1,3) + DATA{c}.com.lpf05(1:end-1,3)*D.*DATA{c}.com.dot(:,3));
    DATA{c}.correction_PD = g*(P*DATA{c}.com.lpf05(1:end-1,3) + D*DATA{c}.com.dot(:,3));

    DATA{c}.belt.speed.eqn_PxD(1) = 0;
    DATA{c}.belt.speed.eqn_PD(1) = 0;
    for i = 2:length(DATA{c}.correction_PxD)
        DATA{c}.belt.speed.eqn_PxD(i) = ((DATA{c}.correction_PxD(i)/fs)*ae(1)*exp(be(1)*DATA{c}.time(i))) + (DATA{c}.belt.speed.eqn_PxD(i-1)*ae(2)*exp(be(2)*DATA{c}.time(i)));
        DATA{c}.belt.speed.eqn_PD(i) = ((DATA{c}.correction_PD(i)/fs)*ae(1)*exp(be(1)*DATA{c}.time(i))) + (DATA{c}.belt.speed.eqn_PD(i-1)*ae(2)*exp(be(2)*DATA{c}.time(i)));
    end
    
%     subplot(4,2,2*(c-1)+1)
    
    figure(1)
    subplot(4,4,c)
    hold on
    
    window = 1:find(DATA{c}.time > 35, 1, 'first');
    plot(DATA{c}.time(window),DATA{c}.belt.speed.raw(window), 'g','linewidth',2);
    plot(DATA{c}.time(window),DATA{c}.com.lpf05(window,3), 'k','linewidth',2);
    plot(DATA{c}.time(window),DATA{c}.com.dot(window,3), 'k:','linewidth',2);
    plot(DATA{c}.time(window),DATA{c}.com.lpf05(window,3).*DATA{c}.com.dot(window,3), 'r-','linewidth',2);
    if c == 13, legend("beltspd", "com","comdot","com*comdot"); end
    if c <= 4, title(condlabel(c)); xlabel("time (s)"); end
    if ismember(c,[1 5 9])
        ylabel("gain = " + num2str(gain([1 5 9] == c))); 
    end
    
    figure(2)
        subplot(4,4,c)
    hold on
    
    window = 1:find(DATA{c}.time > 35, 1, 'first');
    plot(DATA{c}.time(window),DATA{c}.com.lpf05(window,3)/max(abs(DATA{c}.com.lpf05(window,3))), 'k','linewidth',2);
    plot(DATA{c}.time(window),DATA{c}.com.dot(window,3)/max(abs(DATA{c}.com.dot(window,3))), 'k:','linewidth',2);
    
    xxdot = DATA{c}.com.lpf05(window,3).*DATA{c}.com.dot(window,3);
    plot(DATA{c}.time(window),xxdot/max(abs(xxdot)), 'r-','linewidth',2);
    if c == 13, legend("com","comdot","com*comdot"); end
    if c <= 4, title(condlabel(c)); xlabel("time (s)"); end
    if ismember(c,[1 5 9])
        ylabel("normalized" + newline + "gain = " + num2str(gain([1 5 9] == c))); 
    end
    
%     subplot(4,2,2*(c-1)+2)
% 
% %     plot(DATA{c}.time,DATA{c}.rheel(:,3), 'm');
%     hold on
%     plot(DATA{c}.belt.time,DATA{c}.belt.speed.raw, 'b-', 'linewidth',2);
% %     plot(DATA{c}.time,DATA{c}.com.lpf05(:,3), 'k','linewidth',2);
% %     
% %     plot(DATA{c}.belt.time(1:end-1),DATA{c}.belt.speed.dot, 'r');
% %     plot(DATA{c}.time(1:end-1),DATA{c}.com.dot(:,3), 'k:','linewidth',2);
% %     
%     plot(DATA{c}.time(1:end-1),DATA{c}.correction_PxD*normmag, 'g','linewidth',2);
%     plot(DATA{c}.time(1:end-1),DATA{c}.correction_PD*normmag, 'g--','linewidth',2);
% %     plot(DATA{c}.time(1:end-2),DATA{c}.com.ddot(:,3), 'k-.','linewidth',1);
% %     
%     plot(DATA{c}.time(1:end-1),DATA{c}.belt.speed.eqn_PxD, 'r-','linewidth',2);
%     plot(DATA{c}.time(1:end-1),DATA{c}.belt.speed.eqn_PD, 'r--','linewidth',2);
% %     legend("heel", "belt speed", "com", "belt speed dot", "comdot", "correction", "comddot", "belt speed eqn")

end

figure(2)
plot(treadmill.SwayActual, 'k'), hold on
plot(treadmill.RKNEa/100, 'g'), hold on
plot(treadmill.LKNEa/100, 'r')

figure(3)
plot(treadmill.SwayActual, 'k'), hold on
plot(markers.LHEEPosZ, 'g'), hold on
plot(markers.RHEEPosZ, 'r')
ylim([-0.5 0.5])
