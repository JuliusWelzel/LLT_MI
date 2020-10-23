% Get RTs of LLT data
% Clean Overview matrix

PATHIN_RTs = [MAIN '02_data\03_RTs\'];
PATHOUT_RTplot = [MAIN '02_data\03_RTs\group_plot\'];

%Check if PATHOUT-folder is already there; if not: create
if ~isdir(PATHOUT_RTplot)
    mkdir(PATHOUT_RTplot);
end

if exist([PATHIN_RTs 'RT_ALL.mat'],'file') == 2 % update for current evaluation
    load([PATHIN_RTs 'RT_ALL.mat']);
end


%% Delete epochs with invalid speech

for s = 1:length(RT_ALL)
    good_ep = logical(RT_ALL(s).speech_flag);
    RT_ALL(s).SO_ms(~good_ep) = NaN;
end

%% Check Min Max for every participant

for s = 1:length(RT_ALL)
    minRT(s) = min(RT_ALL(s).SO_ms);
    maxRT(s) = max(RT_ALL(s).SO_ms);
end
group = [zeros(1,16),ones(1,23),2*ones(1,22)];

r = figure
plot(1:length(RT_ALL),minRT,1:length(RT_ALL),maxRT)
ylabel 'RTs [ms]'
xlabel 'SUBJ'
ylim ([0 9000])
legend ('minRT','maxRT')
hline ([300,8000],'k:',{'300ms','8000ms'})

%% Plot overview RTs of all participants