%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  RT TIMES IN LLT GROUP PLOT RESULTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Julius Welzel
% University of Oldenburg, 2018

%% Setup
clc;clear all;close all;
MAIN = ['G:\nic_LLT_all\ana03_objectbased\'];
addpath(genpath('G:\nic_LLT_all'));
PATH_RT = [MAIN 'ana03_RTs\ana02_RT_final\'];

load ([PATH_RT 'RT_ALL.mat']); 


%% Sort per group
cs = 1; co = 1; cy = 1;

for i = 1:size(RT_ALL,2)
    if i < 18
        RTs(1).group(cs) = RT_ALL(i);
        cs = cs+1;
    elseif i < 37 && i >17
        RTs(2).group(co) = RT_ALL(i);
        co = co+1;
     elseif i > 36 
        RTs(3).group(cy) = RT_ALL(i);
        cy = cy+1;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       GET ALL RTs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% group: 1 stroke, 2 old, 3 young

% hands
% mh = medial hand, lh = lateral hand
ALL_STROKE = [RTs(1).group.SO_ms];
mh_STROKE = ALL_STROKE([RTs(1).group.con_ind]== 1 | [RTs(1).group.con_ind]== 2 & [RTs(1).group.stim_ind]== 1); % medial hands
lh_STROKE = ALL_STROKE([RTs(1).group.con_ind]== 1 | [RTs(1).group.con_ind]== 2 & [RTs(1).group.stim_ind]== 2); % lateral hands

ALL_OLD = [RTs(2).group.SO_ms];
mh_OLD = ALL_OLD([RTs(2).group.con_ind]== 1 | [RTs(2).group.con_ind]== 2 & [RTs(2).group.stim_ind]== 1); % medial hands
lh_OLD = ALL_OLD([RTs(2).group.con_ind]== 1 | [RTs(2).group.con_ind]== 2 & [RTs(2).group.stim_ind]== 2); % lateral hands

ALL_YOUNG = [RTs(3).group.SO_ms];
mh_YOUNG = ALL_YOUNG([RTs(3).group.con_ind]== 1 | [RTs(3).group.con_ind]== 2 & [RTs(3).group.stim_ind]== 1); % medial hands
lh_YOUNG = ALL_YOUNG([RTs(3).group.con_ind]== 1 | [RTs(3).group.con_ind]== 2 & [RTs(3).group.stim_ind]== 2); % lateral hands


%feet
% mh = medial hand, lh = lateral hand
mf_STROKE = ALL_STROKE([RTs(1).group.con_ind]== 3 | [RTs(1).group.con_ind]== 4 & [RTs(1).group.stim_ind]== 1); % medial hands
lf_STROKE = ALL_STROKE([RTs(1).group.con_ind]== 3 | [RTs(1).group.con_ind]== 4 & [RTs(1).group.stim_ind]== 2); % lateral hands

mf_OLD = ALL_OLD([RTs(2).group.con_ind]== 3 | [RTs(2).group.con_ind]== 4 & [RTs(2).group.stim_ind]== 1); % medial hands
lf_OLD = ALL_OLD([RTs(2).group.con_ind]== 3 | [RTs(2).group.con_ind]== 4 & [RTs(2).group.stim_ind]== 2); % lateral hands

mf_YOUNG = ALL_YOUNG([RTs(3).group.con_ind]== 3 | [RTs(3).group.con_ind]== 4 & [RTs(3).group.stim_ind]== 1); % medial hands
lf_YOUNG = ALL_YOUNG([RTs(3).group.con_ind]== 3 | [RTs(3).group.con_ind]== 4 & [RTs(3).group.stim_ind]== 2); % lateral hands


%% accuracies

acc_STROKE = [RTs(1).group.acc];
acc_mh_STROKE = acc_STROKE([RTs(1).group.con_ind]== 1 | [RTs(1).group.con_ind]== 2 & [RTs(1).group.stim_ind]== 1); % medial hands
acc_lh_STROKE = acc_STROKE([RTs(1).group.con_ind]== 1 | [RTs(1).group.con_ind]== 2 & [RTs(1).group.stim_ind]== 2); % lateral hands
acc_mf_STROKE = acc_STROKE([RTs(1).group.con_ind]== 3 | [RTs(1).group.con_ind]== 4 & [RTs(1).group.stim_ind]== 1); % medial feet
acc_lf_STROKE = acc_STROKE([RTs(1).group.con_ind]== 3 | [RTs(1).group.con_ind]== 4 & [RTs(1).group.stim_ind]== 2); % lateral feet


acc_OLD = [RTs(2).group.acc];
acc_mh_OLD = acc_OLD([RTs(2).group.con_ind]== 1 | [RTs(2).group.con_ind]== 2 & [RTs(2).group.stim_ind]== 1); % medial hands
acc_lh_OLD = acc_OLD([RTs(2).group.con_ind]== 1 | [RTs(2).group.con_ind]== 2 & [RTs(2).group.stim_ind]== 2); % lateral hands
acc_mf_OLD = acc_OLD([RTs(2).group.con_ind]== 3 | [RTs(2).group.con_ind]== 4 & [RTs(2).group.stim_ind]== 1); % medial feet
acc_lf_OLD = acc_OLD([RTs(2).group.con_ind]== 3 | [RTs(2).group.con_ind]== 4 & [RTs(2).group.stim_ind]== 2); % lateral feet

acc_YOUNG = [RTs(3).group.acc];
acc_mh_YOUNG = acc_YOUNG([RTs(3).group.con_ind]== 1 | [RTs(3).group.con_ind]== 2 & [RTs(3).group.stim_ind]== 1); % medial hands
acc_lh_YOUNG = acc_YOUNG([RTs(3).group.con_ind]== 1 | [RTs(3).group.con_ind]== 2 & [RTs(3).group.stim_ind]== 2); % lateral hands
acc_mf_YOUNG = acc_YOUNG([RTs(3).group.con_ind]== 3 | [RTs(3).group.con_ind]== 4 & [RTs(3).group.stim_ind]== 1); % medial feet
acc_lf_YOUNG = acc_YOUNG([RTs(3).group.con_ind]== 3 | [RTs(3).group.con_ind]== 4 & [RTs(3).group.stim_ind]== 2); % lateral feet




for g = 1:3 % group: 1 stroke, 2 old, 3 young
    ACC_hand(g).group = {};
    ACC_feet(g).group = {};
    
    for sub = 1:size(RTs(g).group,2)
%         acc_hand
        mh = RTs(g).group(sub).acc([RTs(g).group(sub).con_ind]== 1 | [RTs(g).group(sub).con_ind]== 2 & [RTs(g).group(sub).stim_ind]== 1);% medial hands
        ACC_hand(g).group(end+1).med = mean(mh);
        lh = RTs(g).group(sub).acc([RTs(g).group(sub).con_ind]== 1 | [RTs(g).group(sub).con_ind]== 2 & [RTs(g).group(sub).stim_ind]== 2);% lateral hands
        ACC_hand(g).group(end+1).lat = mean(lh);
%         acc_feet
        mf = RTs(g).group(sub).acc([RTs(g).group(sub).con_ind]== 3 | [RTs(g).group(sub).con_ind]== 4 & [RTs(g).group(sub).stim_ind]== 1);% medial feet
        ACC_feet(g).group(end+1).med = mean(mf);
        lf = RTs(g).group(sub).acc([RTs(g).group(sub).con_ind]== 3 | [RTs(g).group(sub).con_ind]== 4 & [RTs(g).group(sub).stim_ind]== 2);% lateral feet
        ACC_feet(g).group(end+1).lat = mean(lf);
    end
end

%hand
oh_MedM = mean([ACC_hand(2).group.med]);
oh_MedSD = std([ACC_hand(2).group.med])/sqrt(length([ACC_hand(2).group.med]));
oh_LatM = mean([ACC_hand(2).group.lat]);
oh_LatSD = std([ACC_hand(2).group.lat])/sqrt(length([ACC_hand(2).group.lat]));

yh_MedM = mean([ACC_hand(3).group.med]);
yh_MedSD = std([ACC_hand(3).group.med])/sqrt(length([ACC_hand(3).group.med]));
yh_LatM = mean([ACC_hand(3).group.lat]);
yh_LatSD = std([ACC_hand(3).group.lat])/sqrt(length([ACC_hand(3).group.lat]));

%feet
of_MedM = mean([ACC_feet(2).group.med]);
of_MedSD = std([ACC_feet(2).group.med])/sqrt(length([ACC_feet(2).group.med]));
of_LatM = mean([ACC_feet(2).group.lat]);
of_LatSD = std([ACC_feet(2).group.lat])/sqrt(length([ACC_feet(2).group.lat]));

yf_MedM = mean([ACC_feet(3).group.med]);
yf_MedSD = std([ACC_feet(3).group.med])/sqrt(length([ACC_feet(3).group.med]));
yf_LatM = mean([ACC_feet(3).group.lat]);
yf_LatSD = std([ACC_feet(3).group.lat])/sqrt(length([ACC_feet(3).group.lat]));



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   RAINPLOTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rainplots Hand
hand = figure;
[cb] = cbrewer('qual','Set3',12,'pchip');

subplot(2,3,[1 2])
h1= raincloud_plot('X', lh_OLD, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1,'bandwidth',1 ,'density_type', 'rash');
h2= raincloud_plot('X', lh_YOUNG, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
    'box_col_match', 1,'bandwidth',1, 'density_type', 'rash');
legend([h1{1} h2{1}], {'OLD', 'YOUNG'})
title('Lateral hand RT')
% set(gca,'XLim', [0 40]);
box off
xlabel '[ms]'

subplot(2,3,[4 5])
h3= raincloud_plot('X', mh_OLD, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1,'bandwidth',1 ,'density_type', 'rash');
h4= raincloud_plot('X', mh_YOUNG, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
    'box_col_match', 1,'bandwidth',1, 'density_type', 'rash');
legend([h3{1} h4{1}], {'OLD', 'YOUNG'})
title('Medial hand RT')
% set(gca,'XLim', [0 40]);
box off
xlabel '[ms]'

subplot(2,3,[3 6])

hand_series = [oh_MedM yh_MedM; oh_LatM yh_LatM];
hand_error = [oh_MedSD yh_MedSD; oh_LatSD yh_LatSD];
bh = bar(hand_series);
bh(1).FaceColor = cb(1,:);
bh(2).FaceColor = cb(4,:);
hold on

ngroups = size(hand_series, 1);
nbars = size(hand_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, hand_series(:,i), hand_error(:,i), 'k', 'linestyle', 'none');
end

legend( {'OLD', 'YOUNG'})
xticklabels(gca,{ 'Medial feet', 'Lateral feet'});
ylim ([0.80 0.95]);
title 'Correct %'


%% Rainplots feet
feet = figure;
[cb] = cbrewer('qual','Set3',12,'pchip');

subplot(2,3,[1 2])
h1= raincloud_plot('X', lf_OLD, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1,'bandwidth',1 ,'density_type', 'rash');
h2= raincloud_plot('X', lf_YOUNG, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
    'box_col_match', 1,'bandwidth',1, 'density_type', 'rash');
legend([h1{1} h2{1}], {'OLD', 'YOUNG'})
title('Lateral feet RT')
% set(gca,'XLim', [0 40]);
box off
xlabel '[ms]'

subplot(2,3,[4 5])
h3= raincloud_plot('X', mf_OLD, 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1,'bandwidth',1 ,'density_type', 'rash');
h4= raincloud_plot('X', mf_YOUNG, 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
    'box_col_match', 1,'bandwidth',1, 'density_type', 'rash');
legend([h3{1} h4{1}], {'OLD', 'YOUNG'})
title('Medial feet RT')
% set(gca,'XLim', [0 40]);
box off
xlabel '[ms]'

subplot(2,3,[3 6])

feet_series = [of_MedM yf_MedM; of_LatM yf_LatM];
feet_error = [of_MedSD yf_MedSD; of_LatSD yf_LatSD];
bf = bar(feet_series);
bf(1).FaceColor = cb(1,:);
bf(2).FaceColor = cb(4,:);

hold on

ngroups = size(feet_series, 1);
nbars = size(feet_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, feet_series(:,i), feet_error(:,i), 'k', 'linestyle', 'none');
end

legend( {'OLD', 'YOUNG'})
xticklabels(gca,{ 'Medial feet', 'Lateral feet'});
ylim ([0.80 0.95]);
title 'Correct %'


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   VIOLINPLOTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% lat_med = NaN(2,23,3);
% 
% for g = 1:3 % group: 1 stroke, 2 old, 3 young
%     for sub = 1:size(RTs(g).group,2)
%         med = RTs(g).group(sub).SO_ms(RTs(g).group(sub).con_ind == 1); % medial first row
%         lat_med(1,sub,g) = nanmean(med(11:end));
%         lat = RTs(g).group(sub).SO_ms(RTs(g).group(sub).con_ind == 2); % lateral second row
%         lat_med(2,sub,g) = nanmean(lat(11:end));
%     end
% end
% 
% %% group 
% lSTROKE = lat_med(2,:,1);
% mSTROKE = lat_med(1,:,1);
% 
% lOLD = lat_med(2,:,2);
% mOLD = lat_med(1,:,2);
% 
% lYOUNG = lat_med(2,:,3);
% mYOUNG = lat_med(1,:,3);
% 
% 
% %% plot all RTs
% label = {'lateral','medial'};
% group = {'Stroke','Old','Young'};
% group_idx_med = [ones(1,length(ma_STROKE)),2*ones(1,length(ma_OLD)),3*ones(1,length(ma_YOUNG))];
% group_idx_lat = [ones(1,length(la_STROKE)),2*ones(1,length(la_OLD)),3*ones(1,length(la_YOUNG))];
% 
% figure;
% violinplot([ma_STROKE, ma_OLD, ma_YOUNG]',group_idx_med)
% ylabel '[ms]';
% xticklabels(group);
% title('medial hand');
% 
% figure;
% violinplot([la_STROKE, la_OLD, la_YOUNG]',group_idx_lat)
% xticklabels(group);
% title('lateral hand');
% 


%% Publish

publish('RT_plots.m')









