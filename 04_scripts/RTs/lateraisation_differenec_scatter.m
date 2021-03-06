%% Get RTs to replicate deSimione 2013

clc;clear all;close all;

MAIN = ['E:\nic_LLT_all\ana03_objectbased\ana03_RTs\ana02_RT_final\'];
addpath(genpath('G:\nic_LLT_all\'));
addpath(genpath(MAIN));

load 'RT_all.mat'



[ zero.sixty zero.onet zero.twofor zero.three ] = getRT_con( RT_ALL, [0]);
[ sixty.sixty sixty.onet sixty.twofor sixty.three ] = getRT_con( RT_ALL, [60]);
[ three.sixty three.onet three.twofor three.three ] = getRT_con( RT_ALL, [300]);

[ all.sixty all.onet all.twofor all.three ] = getRT_con( RT_ALL, [0 60 300]);

%% PLOTS for all in depth rotation

sub_c = 1:4;
sub_r = 1:3;
sub_all_c = [4,4,4,4];
% Plot RTs 
figure;

% 0 degree
% scatter
subplot(max(sub_r),max(sub_all_c),[1:sub_c(2)])
sc_o_l = {zero.sixty.lh_o;zero.onet.lh_o;zero.twofor.lh_o;zero.three.lh_o};
sc_o_r = {zero.sixty.rh_o;zero.onet.rh_o;zero.twofor.rh_o;zero.three.rh_o};

for n = 1:4
    scatter(n*ones(length(sc_o_l{n,:}),1),sc_o_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','r', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_o_r{n,:}),1),sc_o_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','g', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_o sem_lh_o] = mean_SEM({zero.sixty.lh_o;zero.onet.lh_o;zero.twofor.lh_o;zero.three.lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','MarkerSize',10,'color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({zero.sixty.rh_o;zero.onet.rh_o;zero.twofor.rh_o;zero.three.rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','MarkerSize',10,'color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 5000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Old 0� '

% scatter
subplot(max(sub_r),max(sub_all_c),[sub_c(3):sub_c(4)])

sc_y_l = {zero.sixty.lh_y;zero.onet.lh_y;zero.twofor.lh_y;zero.three.lh_y};
sc_y_r = {zero.sixty.rh_y;zero.onet.rh_y;zero.twofor.rh_y;zero.three.rh_y};

for n = 1:4
    scatter(n*ones(length(sc_y_l{n,:}),1),sc_y_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','m', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_y_r{n,:}),1),sc_y_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','k', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_y sem_lh_y] = mean_SEM({zero.sixty.lh_y;zero.onet.lh_y;zero.twofor.lh_y;zero.three.lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','MarkerSize',10,'color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({zero.sixty.rh_y;zero.onet.rh_y;zero.twofor.rh_y;zero.three.rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','MarkerSize',10,'color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Young 0�'

% 60 degree
% scatter
subplot(max(sub_r),max(sub_all_c),[sub_c(1):sub_c(2)]+max(sub_all_c))
sc_o_l = {sixty.sixty.lh_o;sixty.onet.lh_o;sixty.twofor.lh_o;sixty.three.lh_o};
sc_o_r = {sixty.sixty.rh_o;sixty.onet.rh_o;sixty.twofor.rh_o;sixty.three.rh_o};

for n = 1:4
    scatter(n*ones(length(sc_o_l{n,:}),1),sc_o_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','r', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_o_r{n,:}),1),sc_o_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','g', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_o sem_lh_o] = mean_SEM({sixty.sixty.lh_o;sixty.onet.lh_o;sixty.twofor.lh_o;sixty.three.lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','MarkerSize',10,'color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({sixty.sixty.rh_o;sixty.onet.rh_o;sixty.twofor.rh_o;sixty.three.rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','MarkerSize',10,'color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 5000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Old 60� '

% scatter
subplot(max(sub_r),max(sub_all_c),[sub_c(3):sub_c(4)]+max(sub_all_c))

sc_y_l = {zero.sixty.lh_y;zero.onet.lh_y;zero.twofor.lh_y;zero.three.lh_y};
sc_y_r = {zero.sixty.rh_y;zero.onet.rh_y;zero.twofor.rh_y;zero.three.rh_y};

for n = 1:4
    scatter(n*ones(length(sc_y_l{n,:}),1),sc_y_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','m', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_y_r{n,:}),1),sc_y_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','k', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_y sem_lh_y] = mean_SEM({sixty.sixty.lh_y;sixty.onet.lh_y;sixty.twofor.lh_y;sixty.three.lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','MarkerSize',10,'color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({sixty.sixty.rh_y;sixty.onet.rh_y;sixty.twofor.rh_y;sixty.three.rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','MarkerSize',10,'color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 5000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Young 60�'

% 300 degree
% scatter 
subplot(max(sub_r),max(sub_all_c),[sub_c(1):sub_c(2)]+2*max(sub_all_c))
sc_o_l = {three.sixty.lh_o;three.onet.lh_o;three.twofor.lh_o;three.three.lh_o};
sc_o_r = {three.sixty.rh_o;three.onet.rh_o;three.twofor.rh_o;three.three.rh_o};

for n = 1:4
    scatter(n*ones(length(sc_o_l{n,:}),1),sc_o_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','r', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_o_r{n,:}),1),sc_o_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','g', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_o sem_lh_o] = mean_SEM({three.sixty.lh_o;three.onet.lh_o;three.twofor.lh_o;three.three.lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','MarkerSize',10,'color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({three.sixty.rh_o;three.onet.rh_o;three.twofor.rh_o;three.three.rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','MarkerSize',10,'color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 5000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Old 300� '

% scatter
subplot(max(sub_r),max(sub_all_c),[sub_c(3):sub_c(4)]+2*max(sub_all_c))
sc_y_l = {three.sixty.lh_y;three.onet.lh_y;three.twofor.lh_y;three.three.lh_y};
sc_y_r = {three.sixty.rh_y;three.onet.rh_y;three.twofor.rh_y;three.three.rh_y};

for n = 1:4
    scatter(n*ones(length(sc_y_l{n,:}),1),sc_y_l{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','m', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
    scatter(n*ones(length(sc_y_r{n,:}),1),sc_y_r{n,:},'jitter','on', 'jitterAmount', 0.05, 'MarkerFaceColor','k', 'MarkerEdgeColor','none','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
    hold on;
end

% meanSEM
[me_lh_y sem_lh_y] = mean_SEM({three.sixty.lh_y;three.onet.lh_y;three.twofor.lh_y;three.three.lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','MarkerSize',10,'color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({three.sixty.rh_y;three.onet.rh_y;three.twofor.rh_y;three.three.rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','MarkerSize',10,'color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 5000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
% legend ('Left Hand','Right Hand')
title 'Young 300�'

%% Set general figure option

% All stimuli
%{
% 0,60,300 degree Average plot for both groups
subplot(max(sub_r),max(sub_all_c),[sub_all_c(1):sub_all_c(2),sub_all_c(1)+max(sub_all_c):sub_all_c(2)+max(sub_all_c),sub_all_c(1)+2*max(sub_all_c):sub_all_c(2)+2*max(sub_all_c)])
[me_lh_o sem_lh_o] = mean_SEM({all.sixty.lh_o;all.onet.lh_o;all.twofor.lh_o;all.three.lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({all.sixty.rh_o;all.onet.rh_o;all.twofor.rh_o;all.three.rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Old all� '

subplot(max(sub_r),max(sub_all_c),[sub_all_c(3):sub_all_c(4),sub_all_c(3)+max(sub_all_c):sub_all_c(4)+max(sub_all_c),sub_all_c(3)+2*max(sub_all_c):sub_all_c(4)+2*max(sub_all_c)])
[me_lh_y sem_lh_y] = mean_SEM({all.sixty.lh_y;all.onet.lh_y;all.twofor.lh_y;all.three.lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({all.sixty.rh_y;all.onet.rh_y;all.twofor.rh_y;all.three.rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Young all�'
%}

% 0 degree
%{
%% Extract 0 degree roation horizontally

c = 1;

for sub = 1:length(RT_ALL)
    for st = 10:length(RT_ALL(sub).stim)
        if strcmp(RT_ALL(sub).stim{st,3},'0')  % only take 0 degreee horizontally
            RT_DS(c).deg = str2num(RT_ALL(sub).stim{st,4}); %get the degree of rotation
            RT_DS(c).ID = str2num(RT_ALL(sub).ID{1}(5:6)); %get the Subject ID
            RT_DS(c).limb = RT_ALL(sub).stim{st,1}; %get limb
            RT_DS(c).mod = RT_ALL(sub).con_ind(st); %get the degree of lateralt = 1 /medial = 2
            RT_DS(c).RT = RT_ALL(sub).SO_ms(st); %get RT 
            RT_DS(c).orient = RT_ALL(sub).stim{st,2}; % get orientation of stimulus back/palm

            c = c+1;
        end
        
    end
end



%% Sort RTs by degree
sixty_lh_o = {};
sixty_lh_y = {};
sixty_rh_o = {};
sixty_rh_y = {};

onet_lh_o = {};
onet_lh_y = {};
onet_rh_o = {};
onet_rh_y = {};

twofor_lh_o = {};
twofor_lh_y = {};
twofor_rh_o = {};
twofor_rh_y = {};

three_lh_o = {};
three_lh_y = {};
three_rh_o = {};
three_rh_y = {};

for r = 1:length(RT_DS)
    if RT_DS(r).deg == 60 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        sixty_lh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 60 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        sixty_lh_y{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 60 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        sixty_rh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 60 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        sixty_rh_y{end+1} = RT_DS(r).RT;

    elseif RT_DS(r).deg == 120 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        onet_lh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 120 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        onet_lh_y{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 120 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        onet_rh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 120 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        onet_rh_y{end+1} = RT_DS(r).RT;

    elseif RT_DS(r).deg == 240 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        twofor_lh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 240 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        twofor_lh_y{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 240 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        twofor_rh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 240 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        twofor_rh_y{end+1} = RT_DS(r).RT;
        
    elseif RT_DS(r).deg == 300 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        three_lh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 300 && strcmp(RT_DS(r).limb,'lh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        three_lh_y{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 300 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 30 && RT_DS(r).ID < 60
        three_rh_o{end+1} = RT_DS(r).RT;
    elseif RT_DS(r).deg == 300 && strcmp(RT_DS(r).limb,'rh') && RT_DS(r).ID >= 60 && RT_DS(r).ID < 100
        three_rh_y{end+1} = RT_DS(r).RT;

    end
        
end
        
sixty_lh_o = [sixty_lh_o{:}];
% sixty_lh_o(isoutlier(sixty_lh_o)) = [];
sixty_lh_y = [sixty_lh_y{:}];
sixty_rh_o = [sixty_rh_o{:}];
sixty_rh_y = [sixty_rh_y{:}];

onet_lh_o = [onet_lh_o{:}];
% onet_lh_o(isoutlier(onet_lh_o)) = [];
onet_lh_y = [onet_lh_y{:}];
onet_rh_o = [onet_rh_o{:}];
onet_rh_y = [onet_rh_y{:}];

twofor_lh_o = [twofor_lh_o{:}];
% twofor_lh_o(isoutlier(twofor_lh_o)) = [];
twofor_lh_y = [twofor_lh_y{:}];
twofor_rh_o = [twofor_rh_o{:}];
twofor_rh_y = [twofor_rh_y{:}];

three_lh_o = [three_lh_o{:}];
% three_lh_o(isoutlier(three_lh_o)) = [];
three_lh_y = [three_lh_y{:}];
three_rh_o = [three_rh_o{:}];
three_rh_y = [three_rh_y{:}];


% Plot RTs 
onet_rh_o(end+1) = mean(onet_rh_o);
onet_lh_y (end) = [];

figure;
subplot(1,2,1)
[me_lh_o sem_lh_o] = mean_SEM([sixty_lh_o;onet_lh_o;twofor_lh_o;three_lh_o]);
errorbar(me_lh_o, sem_lh_o,'-sr','color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM([sixty_rh_o;onet_rh_o;twofor_rh_o;three_rh_o]);
errorbar(me_rh_o, sem_rh_o,'-sg','color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Old'

subplot(1,2,2)
[me_lh_y sem_lh_y] = mean_SEM([sixty_lh_y;onet_lh_y;twofor_lh_y;three_lh_y]);
errorbar(me_lh_y, sem_lh_y,'-sm','color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM([sixty_rh_y;onet_rh_y;twofor_rh_y;three_rh_y]);
errorbar(me_rh_y, sem_rh_y,'-sk','color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Young'






%% average plot for group

old_hand = [sixty_lh_o sixty_rh_o;onet_lh_o onet_rh_o;twofor_lh_o twofor_rh_o;three_lh_o three_rh_o];
young_hand = [sixty_lh_y sixty_rh_y;onet_lh_y onet_rh_y;twofor_lh_y twofor_rh_y;three_lh_y three_rh_y];

[m_o_h sem_o_h] = mean_SEM(old_hand);
[m_y_h sem_y_h] = mean_SEM(young_hand);

figure
errorbar(m_o_h, sem_o_h,'-sr','color','r')
hold on;
errorbar(m_y_h, sem_y_h,'-sb','color','b')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Old','Young')


%}

% 60 degree
%{
%% Extract 60 degree roation horizontally

c = 1;

for sub = 1:length(RT_ALL)
    for st = 10:length(RT_ALL(sub).stim)
        if strcmp(RT_ALL(sub).stim{st,3},'60')  % only take 0 degreee horizontally
            RT_DS_six(c).deg = str2num(RT_ALL(sub).stim{st,4}); %get the degree of rotation
            RT_DS_six(c).ID = str2num(RT_ALL(sub).ID{1}(5:6)); %get the Subject ID
            RT_DS_six(c).limb = RT_ALL(sub).stim{st,1}; %get limb
            RT_DS_six(c).mod = RT_ALL(sub).con_ind(st); %get the degree of lateralt = 1 /medial = 2
            RT_DS_six(c).RT = RT_ALL(sub).SO_ms(st); %get RT 
            RT_DS_six(c).orient = RT_ALL(sub).stim{st,2}; % get orientation of stimulus back/palm

            c = c+1;
        end
        
    end
end



%% Sort RTs by degree
sixty_lh_o = {};
sixty_lh_y = {};
sixty_rh_o = {};
sixty_rh_y = {};

onet_lh_o = {};
onet_lh_y = {};
onet_rh_o = {};
onet_rh_y = {};

twofor_lh_o = {};
twofor_lh_y = {};
twofor_rh_o = {};
twofor_rh_y = {};

three_lh_o = {};
three_lh_y = {};
three_rh_o = {};
three_rh_y = {};

for r = 1:length(RT_DS_six)
    if RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        sixty_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        sixty_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        sixty_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        sixty_rh_y{end+1} = RT_DS_six(r).RT;

    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        onet_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        onet_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        onet_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        onet_rh_y{end+1} = RT_DS_six(r).RT;

    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        twofor_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        twofor_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        twofor_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        twofor_rh_y{end+1} = RT_DS_six(r).RT;
        
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        three_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        three_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        three_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        three_rh_y{end+1} = RT_DS_six(r).RT;

    end
        
end
        
sixty_lh_o = [sixty_lh_o{:}];
% sixty_lh_o(isoutlier(sixty_lh_o)) = [];
sixty_lh_y = [sixty_lh_y{:}];
sixty_rh_o = [sixty_rh_o{:}];
sixty_rh_y = [sixty_rh_y{:}];

onet_lh_o = [onet_lh_o{:}];
% onet_lh_o(isoutlier(onet_lh_o)) = [];
onet_lh_y = [onet_lh_y{:}];
onet_rh_o = [onet_rh_o{:}];
onet_rh_y = [onet_rh_y{:}];

twofor_lh_o = [twofor_lh_o{:}];
% twofor_lh_o(isoutlier(twofor_lh_o)) = [];
twofor_lh_y = [twofor_lh_y{:}];
twofor_rh_o = [twofor_rh_o{:}];
twofor_rh_y = [twofor_rh_y{:}];

three_lh_o = [three_lh_o{:}];
% three_lh_o(isoutlier(three_lh_o)) = [];
three_lh_y = [three_lh_y{:}];
three_rh_o = [three_rh_o{:}];
three_rh_y = [three_rh_y{:}];


% Plot RTs 

figure;
subplot(1,2,1)
[me_lh_o sem_lh_o] = mean_SEM({sixty_lh_o;onet_lh_o;twofor_lh_o;three_lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({sixty_rh_o;onet_rh_o;twofor_rh_o;three_rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Old'

subplot(1,2,2)
[me_lh_y sem_lh_y] = mean_SEM({sixty_lh_y;onet_lh_y;twofor_lh_y;three_lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({sixty_rh_y;onet_rh_y;twofor_rh_y;three_rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Young'


% average plot for group

old_hand = {[sixty_lh_o sixty_rh_o];[onet_lh_o onet_rh_o];[twofor_lh_o twofor_rh_o];[three_lh_o three_rh_o]};
young_hand = {[sixty_lh_y sixty_rh_y];[onet_lh_y onet_rh_y];[twofor_lh_y twofor_rh_y];[three_lh_y three_rh_y]};

[m_o_h sem_o_h] = mean_SEM(old_hand);
[m_y_h sem_y_h] = mean_SEM(young_hand);

figure
errorbar(m_o_h, sem_o_h,'-sr','color','r')
hold on;
errorbar(m_y_h, sem_y_h,'-sb','color','b')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Old','Young')

%}

% 300 degree
%{
%% Extract 300 degree roation horizontally

c = 1;

for sub = 1:length(RT_ALL)
    for st = 10:length(RT_ALL(sub).stim)
        if strcmp(RT_ALL(sub).stim{st,3},'300')  % only take 0 degreee horizontally
            RT_DS_six(c).deg = str2num(RT_ALL(sub).stim{st,4}); %get the degree of rotation
            RT_DS_six(c).ID = str2num(RT_ALL(sub).ID{1}(5:6)); %get the Subject ID
            RT_DS_six(c).limb = RT_ALL(sub).stim{st,1}; %get limb
            RT_DS_six(c).mod = RT_ALL(sub).con_ind(st); %get the degree of lateralt = 1 /medial = 2
            RT_DS_six(c).RT = RT_ALL(sub).SO_ms(st); %get RT 
            RT_DS_six(c).orient = RT_ALL(sub).stim{st,2}; % get orientation of stimulus back/palm

            c = c+1;
        end
        
    end
end



%% Sort RTs by degree
sixty_lh_o = {};
sixty_lh_y = {};
sixty_rh_o = {};
sixty_rh_y = {};

onet_lh_o = {};
onet_lh_y = {};
onet_rh_o = {};
onet_rh_y = {};

twofor_lh_o = {};
twofor_lh_y = {};
twofor_rh_o = {};
twofor_rh_y = {};

three_lh_o = {};
three_lh_y = {};
three_rh_o = {};
three_rh_y = {};

for r = 1:length(RT_DS_six)
    if RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        sixty_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        sixty_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        sixty_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 60 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        sixty_rh_y{end+1} = RT_DS_six(r).RT;

    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        onet_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        onet_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        onet_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 120 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        onet_rh_y{end+1} = RT_DS_six(r).RT;

    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        twofor_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        twofor_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        twofor_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 240 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        twofor_rh_y{end+1} = RT_DS_six(r).RT;
        
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        three_lh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'lh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        three_lh_y{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 30 && RT_DS_six(r).ID < 60
        three_rh_o{end+1} = RT_DS_six(r).RT;
    elseif RT_DS_six(r).deg == 300 && strcmp(RT_DS_six(r).limb,'rh') && RT_DS_six(r).ID >= 60 && RT_DS_six(r).ID < 100
        three_rh_y{end+1} = RT_DS_six(r).RT;

    end
        
end
        
sixty_lh_o = [sixty_lh_o{:}];
% sixty_lh_o(isoutlier(sixty_lh_o)) = [];
sixty_lh_y = [sixty_lh_y{:}];
sixty_rh_o = [sixty_rh_o{:}];
sixty_rh_y = [sixty_rh_y{:}];

onet_lh_o = [onet_lh_o{:}];
% onet_lh_o(isoutlier(onet_lh_o)) = [];
onet_lh_y = [onet_lh_y{:}];
onet_rh_o = [onet_rh_o{:}];
onet_rh_y = [onet_rh_y{:}];

twofor_lh_o = [twofor_lh_o{:}];
% twofor_lh_o(isoutlier(twofor_lh_o)) = [];
twofor_lh_y = [twofor_lh_y{:}];
twofor_rh_o = [twofor_rh_o{:}];
twofor_rh_y = [twofor_rh_y{:}];

three_lh_o = [three_lh_o{:}];
% three_lh_o(isoutlier(three_lh_o)) = [];
three_lh_y = [three_lh_y{:}];
three_rh_o = [three_rh_o{:}];
three_rh_y = [three_rh_y{:}];


% Plot RTs 

figure;
subplot(1,2,1)
[me_lh_o sem_lh_o] = mean_SEM({sixty_lh_o;onet_lh_o;twofor_lh_o;three_lh_o});
errorbar(me_lh_o, sem_lh_o,'-sr','color','r')
hold on;
[me_rh_o sem_rh_o] = mean_SEM({sixty_rh_o;onet_rh_o;twofor_rh_o;three_rh_o});
errorbar(me_rh_o, sem_rh_o,'-sg','color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Old'

subplot(1,2,2)
[me_lh_y sem_lh_y] = mean_SEM({sixty_lh_y;onet_lh_y;twofor_lh_y;three_lh_y});
errorbar(me_lh_y, sem_lh_y,'-sm','color','m')
hold on;
[me_rh_y sem_rh_y] = mean_SEM({sixty_rh_y;onet_rh_y;twofor_rh_y;three_rh_y});
errorbar(me_rh_y, sem_rh_y,'-sk','color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Young'


% average plot for group

old_hand = {[sixty_lh_o sixty_rh_o];[onet_lh_o onet_rh_o];[twofor_lh_o twofor_rh_o];[three_lh_o three_rh_o]};
young_hand = {[sixty_lh_y sixty_rh_y];[onet_lh_y onet_rh_y];[twofor_lh_y twofor_rh_y];[three_lh_y three_rh_y]};

[m_o_h sem_o_h] = mean_SEM(old_hand);
[m_y_h sem_y_h] = mean_SEM(young_hand);

figure
errorbar(m_o_h, sem_o_h,'-sr','color','r')
hold on;
errorbar(m_y_h, sem_y_h,'-sb','color','b')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Old','Young')

%}

% All stimuli
%{
%% Extract 0,60 & 300 degree roation horizontally

c = 1;

for sub = 1:length(RT_ALL)
    for st = 10:length(RT_ALL(sub).stim)
            RT_DS_all(c).deg = str2num(RT_ALL(sub).stim{st,4}); %get the degree of rotation
            RT_DS_all(c).ID = str2num(RT_ALL(sub).ID{1}(5:6)); %get the Subject ID
            RT_DS_all(c).limb = RT_ALL(sub).stim{st,1}; %get limb
            RT_DS_all(c).mod = RT_ALL(sub).con_ind(st); %get the degree of lateralt = 1 /medial = 2
            RT_DS_all(c).RT = RT_ALL(sub).SO_ms(st); %get RT 
            RT_DS_all(c).orient = RT_ALL(sub).stim{st,2}; % get orientation of stimulus back/palm

            c = c+1;
        
    end
end


%% Sort RTs by degree
sixty_lh_o_all = {};
sixty_lh_y_all = {};
sixty_rh_o_all = {};
sixty_rh_y_all = {};

onet_lh_o_all = {};
onet_lh_y_all = {};
onet_rh_o_all = {};
onet_rh_y_all = {};

twofor_lh_o_all = {};
twofor_lh_y_all = {};
twofor_rh_o_all = {};
twofor_rh_y_all = {};

three_lh_o_all = {};
three_lh_y_all = {};
three_rh_o_all = {};
three_rh_y_all = {};

for r = 1:length(RT_DS_all)
    if RT_DS_all(r).deg == 60 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        sixty_lh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 60 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        sixty_lh_y_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 60 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        sixty_rh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 60 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        sixty_rh_y_all{end+1} = RT_DS_all(r).RT;

    elseif RT_DS_all(r).deg == 120 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        onet_lh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 120 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        onet_lh_y_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 120 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        onet_rh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 120 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        onet_rh_y_all{end+1} = RT_DS_all(r).RT;

    elseif RT_DS_all(r).deg == 240 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        twofor_lh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 240 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        twofor_lh_y_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 240 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        twofor_rh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 240 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        twofor_rh_y_all{end+1} = RT_DS_all(r).RT;
        
    elseif RT_DS_all(r).deg == 300 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        three_lh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 300 && strcmp(RT_DS_all(r).limb,'lh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        three_lh_y_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 300 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 30 && RT_DS_all(r).ID < 60
        three_rh_o_all{end+1} = RT_DS_all(r).RT;
    elseif RT_DS_all(r).deg == 300 && strcmp(RT_DS_all(r).limb,'rh') && RT_DS_all(r).ID >= 60 && RT_DS_all(r).ID < 100
        three_rh_y_all{end+1} = RT_DS_all(r).RT;

    end
        
end
        
sixty_lh_o_all = [sixty_lh_o_all{:}];
% sixty_lh_o(isoutlier(sixty_lh_o)) = [];
sixty_lh_y_all = [sixty_lh_y_all{:}];
sixty_rh_o_all = [sixty_rh_o_all{:}];
sixty_rh_y_all = [sixty_rh_y_all{:}];

onet_lh_o_all = [onet_lh_o_all{:}];
% onet_lh_o(isoutlier(onet_lh_o)) = [];
onet_lh_y_all = [onet_lh_y_all{:}];
onet_rh_o_all = [onet_rh_o_all{:}];
onet_rh_y_all = [onet_rh_y_all{:}];

twofor_lh_o_all = [twofor_lh_o_all{:}];
% twofor_lh_o(isoutlier(twofor_lh_o)) = [];
twofor_lh_y_all = [twofor_lh_y_all{:}];
twofor_rh_o_all = [twofor_rh_o_all{:}];
twofor_rh_y_all = [twofor_rh_y_all{:}];

three_lh_o_all = [three_lh_o_all{:}];
% three_lh_o(isoutlier(three_lh_o)) = [];
three_lh_y_all = [three_lh_y_all{:}];
three_rh_o_all = [three_rh_o_all{:}];
three_rh_y_all = [three_rh_y_all{:}];

%% Plot RTs 

figure;
subplot(1,2,1)
[me_lh_o_all sem_lh_o_all] = mean_SEM({sixty_lh_o_all;onet_lh_o_all;twofor_lh_o_all;three_lh_o_all});
errorbar(me_lh_o_all, sem_lh_o_all,'-sr','color','r')
hold on;
[me_rh_o_all sem_rh_o_all] = mean_SEM({sixty_rh_o_all;onet_rh_o_all;twofor_rh_o_all;three_rh_o_all});
errorbar(me_rh_o_all, sem_rh_o_all,'-sg','color','g')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Old'

subplot(1,2,2)
[me_lh_y_all sem_lh_y_all] = mean_SEM({sixty_lh_y_all;onet_lh_y_all;twofor_lh_y_all;three_lh_y_all});
errorbar(me_lh_y_all, sem_lh_y_all,'-sm','color','m')
hold on;
[me_rh_y_all sem_rh_y_all] = mean_SEM({sixty_rh_y_all;onet_rh_y_all;twofor_rh_y_all;three_rh_y_all});
errorbar(me_rh_y_all, sem_rh_y_all,'-sk','color','k')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Left Hand','Right Hand')
title 'Young'






%% average plot for group

old_hand = [sixty_lh_o sixty_rh_o;onet_lh_o onet_rh_o;twofor_lh_o twofor_rh_o;three_lh_o three_rh_o];
young_hand = [sixty_lh_y sixty_rh_y;onet_lh_y onet_rh_y;twofor_lh_y twofor_rh_y;three_lh_y three_rh_y];

[m_o_h sem_o_h] = mean_SEM(old_hand);
[m_y_h sem_y_h] = mean_SEM(young_hand);

figure
errorbar(m_o_h, sem_o_h,'-sr','color','r')
hold on;
errorbar(m_y_h, sem_y_h,'-sb','color','b')

xlim([0 5])
xlabel 'Degree in rotation'
ylim([1000 2000])
ylabel 'RT[ms]'
xticks([1:4])
xticklabels({'60�','120�','240�','300�'})
legend ('Old','Young')

%}
