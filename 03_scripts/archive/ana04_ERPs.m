%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       ERP of LLT stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ERP analysis of LLT group data 
% Data: nic_LLT (Niclas Braun, University of Oldenburg)
% Author: Julius Welzel, 0Mareike Daeglau & Catharina Zich 
% Supervisor: Cornelia Kranczioch (University of Oldenburg)


global MAIN
PATHIN_ERP = [MAIN '02_data\04_ERPs\'];

PATHOUT_ERP = [MAIN '02_data\04_ERPs\'];
PATHOUT_ERPplot = [PATHOUT_ERP 'ERP_plots\'];


%Check if PATHOUT-folder is already there; if not: create

if ~isdir(PATHOUT_ERP)
    mkdir(PATHOUT_ERP);
elseif ~isdir(PATHOUT_ERPplot)
    mkdir(PATHOUT_ERPplot)
end

% document potential problems with a subject
docError = {};   

% load RT data for check for valid trials
if exist([MAIN '02_data\03_RTs\RT_ALL.mat'],'file') == 2 % update for current evaluation
    load([MAIN '02_data\03_RTs\RT_ALL.mat']);
end

%% Gather all available datasets with clean EEG data
list = dir(fullfile([PATHIN_ERP '*ep_filt.set']));
SUBJ = str2double(extractBetween({list.name},'SUBJ','_'));

%%
for sub = SUBJ(old|young)
        EEG = pop_loadset('filename',['SUBJ' num2str(sub) '_ep_filt.set'],'filepath',PATHIN_ERP);
        ERP(sub).ID = EEG.ID;
        ERP(sub).data = EEG.data;
end


%%
old = SUBJ<70 & SUBJ>=30;
young = SUBJ>= 70;

LH_O_P = {};
con = {'LH','RH','LF','RF'};


for c = 1:length(con)
    % get relevant trials in Old group for all 4 conditions
    eval ([con{c} '_O = {}']);
    for s = SUBJ(old)
        s_rt = find(strcmp({RT_ALL.ID},[ERP(s).ID]));
        trials = RT_ALL(s_rt).con_ind(11:end)==1 & RT_ALL(s_rt).acc(11:end)==1 & RT_ALL(s_rt).speech_flag(11:end)==1;
        eval([con{c} '_O{end+1} = [ERP(s).data(:,:,trials)];']);
    end
    
    % put trials back in 3D matrix
    eval([con{c} '_O = cat(3,' con{c} '_O{:});'])

end



for c = 1:length(con)
    % get relevant trials in Old group for all 4 conditions
    eval ([con{c} '_Y = {}']);
    for s = SUBJ(young)
        s_rt = find(strcmp({RT_ALL.ID},[ERP(s).ID]));
        trials = RT_ALL(s_rt).con_ind(11:end)==1 & RT_ALL(s_rt).acc(11:end)==1 & RT_ALL(s_rt).speech_flag(11:end)==1;
        eval([con{c} '_Y{end+1} = [ERP(s).data(:,:,trials)];']);
    end
    
    % put trials back in 3D matrix
    eval([con{c} '_Y = cat(3,' con{c} '_Y{:});'])

end

%% Plot ERP
e = figure;
subplot(2,1,1)
plot(EEG.times,mean(LH_O,3),'color',cc(1,:))


%%%%%%%%%%%%%%%%%%%%% Graphics Settings - can be customized %%%%%%%%%%%%%%%%%%
%
LINEWIDTH     = 1;     % data line widths (can be non-integer)
FONTSIZE      = 10;      % font size to use for labels
CHANFONTSIZE  = 7;       % font size to use for channel names
TICKFONTSIZE  = 8;       % font size to use for axis labels
TITLEFONTSIZE = 12;      % font size to use for the plot title
PLOT_WIDTH    = 0.95;     % 0.75, width and height of plot array on figure
PLOT_HEIGHT   = 0.88;    % 0.88
gcapos = get(gca,'Position'); axis off;
PLOT_WIDTH    = gcapos(3)*PLOT_WIDTH; % width and height of gca plot array on gca
PLOT_HEIGHT   = gcapos(4)*PLOT_HEIGHT;
curfig = gcf;            % learn the current graphic figure number

 % read the channel location file
 % ------------------------------
if isstruct(EEG.chanlocs)
    nonemptychans = cellfun('isempty', { EEG.chanlocs.theta });
    nonemptychans = find(~nonemptychans);
    [tmp EEG.channames Th Rd] = readlocs(EEG.chanlocs(nonemptychans));
    EEG.channames = strvcat({ EEG.chanlocs.labels });
else
    [tmp EEG.channames Th Rd] = readlocs(EEG.chanlocs);
    EEG.channames = strvcat(EEG.channames);
    nonemptychans = [1:length(EEG.channames)];
end;
Th = pi/180*Th;                 % convert degrees to radians
Rd = Rd; 

[yvalstmp,xvalstmp] = pol2cart(Th,Rd); % translate from polar to cart. coordinates
xvals(nonemptychans) = xvalstmp;
yvals(nonemptychans) = yvalstmp;

% find position for other channels
% --------------------------------
[EEG.chans,framestotal]=size(EEG.data);           % data size

totalchans = length(EEG.chanlocs);
emptychans = setdiff_bc(1:totalchans, nonemptychans);
totalchans = floor(sqrt(totalchans))+1;
for index = 1:length(emptychans)
    xvals(emptychans(index)) = 0.7+0.2*floor((index-1)/totalchans);
    yvals(emptychans(index)) = -0.4+mod(index-1,totalchans)/totalchans;
end;
EEG.channames = EEG.channames(EEG.chans,:);
xvals     = xvals(EEG.chans);
yvals     = yvals(EEG.chans);


if length(xvals) > 1
    if length(unique(xvals)) > 1
        xvals = (xvals-mean([max(xvals) min(xvals)]))/(max(xvals)-min(xvals)); % recenter
        xvals = gcapos(1)+gcapos(3)/2+PLOT_WIDTH*xvals;   % controls width of plot 
                                                          % array on current axes
    end;
end;
yvals = gcapos(2)+gcapos(4)/2+PLOT_HEIGHT*yvals;  % controls height of plot 
                                                      % array on current axes

