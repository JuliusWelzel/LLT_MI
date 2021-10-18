%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Check RTs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RT analysis of group data
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


PATHIN_RTs  = [MAIN '02_data\03_RTs\'];

% load all ERP data
load ([PATHIN_RTs 'RT_all.mat']);

SUBJ = str2double(extractAfter({RT_ALL.ID},'SUBJ'));

% define groups
idx_old     = SUBJ <70 & SUBJ >= 30;
idx_young   = SUBJ >= 70;
idx_stroke  = SUBJ < 30;
idx_group   = SUBJ;
idx_group(idx_stroke)   = 0;
idx_group(idx_old)      = 1;
idx_group(idx_young)    = 2;

% groups vars for plotting
nms_group   = {'stroke','old','young'}; 

%% plot per group per degree


for s = 1:numel(RT_ALL)
    
    % left hand(first row) all lateral rotations
    rt_sort(1,1,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'lh') & strcmp({RT_ALL(s).stim{:,4}},'60')),'omitnan'); 
    rt_sort(1,2,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'lh') & strcmp({RT_ALL(s).stim{:,4}},'120')),'omitnan'); 
    rt_sort(1,3,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'lh') & strcmp({RT_ALL(s).stim{:,4}},'240')),'omitnan'); 
    rt_sort(1,4,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'lh') & strcmp({RT_ALL(s).stim{:,4}},'300')),'omitnan'); 

    % left hand(second row) all lateral rotations
    rt_sort(2,1,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'rh') & strcmp({RT_ALL(s).stim{:,4}},'60')),'omitnan'); 
    rt_sort(2,2,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'rh') & strcmp({RT_ALL(s).stim{:,4}},'120')),'omitnan'); 
    rt_sort(2,3,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'rh') & strcmp({RT_ALL(s).stim{:,4}},'240')),'omitnan'); 
    rt_sort(2,4,s)  = mean(RT_ALL(s).SO_ms(RT_ALL(s).acc == 1 & strcmp({RT_ALL(s).stim{:,1}},'rh') & strcmp({RT_ALL(s).stim{:,4}},'300')),'omitnan'); 

end

%% plot data

% transform struct to table for better handling
figure

% young
subplot(2,2,1)
ylh = errorbar(mean(rt_sort(1,:,idx_young),3),std(rt_sort(1,:,idx_young),0,3)/sqrt(size(rt_sort(1,:,idx_young),3)),'-ok')
    ylh.MarkerFaceColor = [.3 .3 .3];
    ylh.MarkerSize = 10;
hold on
yrh = errorbar(mean(rt_sort(2,:,idx_young),3),std(rt_sort(2,:,idx_young),0,3)/sqrt(size(rt_sort(2,:,idx_young),3)),'-ok')
    yrh.MarkerFaceColor = [.7 .7 .7];
    yrh.MarkerSize = 10;

    l = legend ({'left hand','right hand'})
    title (l,'Young')
    ylabel 'Time [ms]'
    xticklabels ({'60°','120°','240°','300°'})
    
% old
subplot(2,2,2)
olh = errorbar(mean(rt_sort(1,:,idx_old),3),std(rt_sort(1,:,idx_old),0,3)/sqrt(size(rt_sort(1,:,idx_old),3)),'-ok')
    olh.MarkerFaceColor = [.3 .3 .3];
    olh.MarkerSize = 10;
hold on
orh = errorbar(mean(rt_sort(2,:,idx_old),3),std(rt_sort(2,:,idx_old),0,3)/sqrt(size(rt_sort(2,:,idx_old),3)),'-ok')
    orh.MarkerFaceColor = [.7 .7 .7];
    orh.MarkerSize = 10;

    l = legend ({'left hand','right hand'})
    title (l,'Old')
    ylabel 'Time [ms]'
    xticklabels ({'60°','120°','240°','300°'})
    
subplot(2,2,[3 4])
slh = errorbar(mean(rt_sort(1,:,idx_young),3),std(rt_sort(1,:,idx_young),0,3)/sqrt(size(rt_sort(1,:,idx_young),3)),'-ok')
    slh.MarkerFaceColor = [.3 .3 .3];
    slh.MarkerSize = 10;
hold on
srh = errorbar(mean(rt_sort(2,:,idx_young),3),std(rt_sort(2,:,idx_young),0,3)/sqrt(size(rt_sort(2,:,idx_young),3)),'-ok')
    srh.MarkerFaceColor = [.7 .7 .7];
    srh.MarkerSize = 10;

    l = legend ({'left hand','right hand'})
    title (l,'Stroke')
    ylabel 'Time [ms]'
    xticks ([1:4])
    xticklabels ({'60°','120°','240°','300°'})
