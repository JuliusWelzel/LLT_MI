%% Matlab Befehle
pwd     % identify current path
dbquit  % end debug modus
format longG    % don`t show confusing e^ notation

%return  % return to invoking function
reshape([1 2 3; 4 5 5; 7 8 9], 1,9); % transforme matrix in 1 x 9 array
numel(A) % Number of elements in array or subscripted array expression

%% Plot a selected matrix as an image use (e.g. an erp_image)
figure; imagesc(squeeze(mean(EEG.data,1)))

%% Modulo operator for only choosing every n-th trial:
mod(9,3)==0

% Draw vertical line as at specific sample time point
line([t(i) t(i)],[0 2],'color','r');  % draw a vertical line at time point t

%% Use classify of Matlab
% testFeature1: are all the test feature values of the 1st feature (regardless of class belonging)
% --> each row represents one feature.
% testFeature2: are all the test feature values of the 2nd feature (regardless of class belonging)
% trainFeature1: are all the train feature values of the 1st feature (regardless of class belonging)
% trainFeature2: are all the train feature values of the 2nd feature (regardless of class belonging)
% labelTrian: vector of labels (must have the same row size as testFeature1,2)

[C,err,P,logp,coeff] = classify([testFeature1 testFeature2],[trainFeature1 trainFeature2],...
                                labelTrain,'Quadratic');
   

%% Finding the index value corresponding to a value closest to 0 in an array
%Do you have the Stats toolbox? if you do then do a nearest neighbor search as follows:
location = knnsearch(Temp',30);

%If you don't have stats toolbox then use delauny to do nn search:
tri = delaunayn(Temp');
dsearchn(Temp',tri,30)

%% find-largest-positive-value-next-to-zero/
http://blogs.mathworks.com/loren/2009/07/01/find-largest-positive-value-next-to-zero/
A =  [1, 5, 3, 0, 2, 7, 0, 8, 9, 1 0];
fs = @(A) max(A(imdilate(A == 0, [1 1 1])));
fs(A)

%% Calculate strength of effect
% Als Schätzer für gleiche Gruppengrößen und unterschiedliche Varianzen wurde von Cohen

    d = \frac{\bar{x}_1-\bar{x}_2}{\sqrt{(s_1^2+s_2^2) /2}} 

%% Concatenate a matrix with itself N times
a=(1:4)';
b=repmat(a,1,7)
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/251653

%% What is a covariance matrix and how to calculate it
% each row is one observation, each column is one variable!
x = [   1     1     4;
        2     2    4;
        1     3    4;
        2     4    4];
cov(x)

% gives
        var 1     var 2          var 3
var1    0.3333    0.3333         0
var2    0.3333    1.6667         0
var3    0         0              0


% matrix diagnoal gives variance of each variable within itself.
% values above and below diagonal are mirrored (redundant)
% These values give the covariance (relation) of each variable pair
% --> if > 0 positive relation
% --> if < 0 negative relation
% --> if == 0 no relation.
% --> Stength of relation cannot be inferred.



%% Sort an array and get the indices:
[x,i] = sort([10:-1:5])

%% Get figure handles of all open figures
 figHandles = get(0,'Children');
 close(figHandles(2:end)) % close all figures apart from the first one

%% Calculate multiple ttests (e.g. for plotting significance of signal difference)
% each column represents one distribution (one t-test)
a =[1 1 2    3 1; 
    2 1 2    3 2; 
    3 1 2.5  3 3; 
    4 1 2    4 3; 
    5 1 3.5  5 5];

b = [1 1 1 1 1; 
    2 2 2 2 2; 
    3 3 3 3 3; 
    4 4 4 4 4; 
    5 5 5 5 5];


[h,p] = ttest(a,b)

%% generates linearly spaced vectors. 
% The linspace function generates linearly spaced vectors. 
% It is similar to the colon operator ":" but gives direct control over the number of points and always includes the endpoints.
y = linspace(x1,x2) % returns a row vector with 100 linearly 
                    % spaced points in the interval [x1,x2].
                    
y = linspace(x1,x2,N) % returns N linearly spaced points.

%% more precise timing:
% use  java.lang.Thread.sleep(1); instead of pause() command.
% http://undocumentedmatlab.com/blog/pause-for-the-better

%% interpolate nan
a= [1 2 nan 3 4 5 nan nan 6 7 8 9]
inpaint_nans(a)

%% resample time series
y = resample(x,1,2); % resample to half the sample rate
y = resample(x,2,1); % resample to double the sample rate

%% fast refreshing of a real time block:
% http://stackoverflow.com/questions/13102654/how-should-i-update-the-data-of-a-plot-in-matlab

%% make bar plot with errors:
values = randn(4,3);         % random y values (4 groups of 3 parameters)
errY = zeros(4,3);
errY(:,:,1) = 0.1.*values;   % 10% lower error

figure();
    h = barwitherr(values, errY);% Plot with errorbars
 
    set(gca,'XTickLabel',{'Group 1','Group 2', 'Group 3','Group 4'},'fontsize',13)
    legend({'all', '100', '400', '700'},'fontsize',11);
    ylabel('Estimation error (Percentage)');
    set(h(1),'FaceColor','k');
    %title('Intentional Binding ([actual-estimate]/actual)','fontsize',14);
    title('Intentional Binding actual-estimate','fontsize',15);
    set(gcf, 'Units','centimeters','Position', [1, 1, 20, 13]);

%% Gentle way interrupting while with button press:
k=0;
while ~k
    k=waitforbuttonpress;
    if ~strcmp(get(gcf,'currentcharacter'),'5');
        k=0;
    end
end

%% Counter number of occurenes of pairs that occur in two arrays:
A = [1 6 1 10 1 7 1 9  3 6 0 0 -2 -2];
B = [7 2 3  5 6 8 7 9 10 2 0 0 1 0];

%// Find the unique pairs
AB  = [A;B].';
ABu = unique(AB, 'rows');

%// Count the number of occurrences of each unique pair
%// (pretty sure there's a better way to do this...)     
C = arrayfun(@(ii) ...
    [ABu(ii,:) sum(all(bsxfun(@eq, AB, ABu(ii,:)),2))], 1:size(ABu,1), ...
    'UniformOutput', false);
C = cat(1,C{:});


    
%% How to color bar graphs individually
% Color
data = 1:12
clrCerulean = [0.0, 0.48, 0.65];
clrOrangeRed = [1.0, 0.27, 0.0];
clrOliveGreen = [0.33, 0.42, 0.18];
 
% Plot bar graph and color it
% Bar graph figure
hBar = bar(data);
 
% Get a handle to the children
hBarChildren = get(hBar, 'Children');
 
% Set the colors we want to use
myBarColors = [clrCerulean; clrOrangeRed; clrOliveGreen];
 
% This defines which bar will be using which index of "myBarColors", i.e. the first
% two bars will be colored in "clrCerulean", the next 6 will be colored in "clrOrangeRed"
% and the last 4 bars will be colored in "clrOliveGreen"
index = [1 1 2 2 2 2 2 2 3 3 3 3];
 
% Set the index as CData to the children
set(hBarChildren, 'CData', index);
 
% And set the custom colormap. Takes care of everything else
colormap(myBarColors);

%% add a supertitle in subplots:
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Distributions of Intentional Binding', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%% subtract a vector from a matrix rowwise:
a = [1 2 3; 4 5 6; 7 8 9]; % my matrix
b = [1 2 3];                % my vector
c = bsxfun(@minus,a,b)      % my resulting matrix;

%% make a scotter plot (correlation between two variables)
scatter(x,y,[],'b'); % plot scatterplot in blue
lsline; % best fit line


%% install and open eeglab
cp = 'C:\Users\nib\Desktop\Eyetracking';
addpath([cp '/eeglab10.2.2.4b/']);
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% load and save files
EEG = pop_biosig('C:\Users\nib\Desktop\Eyetracking\stefanSync.gdf');
EEG = pop_loadset( filename, filepath);
EEG = pop_saveset(EEG, 'filename',['stefanSync.set'], 'filepath',cp);
eeglab redraw;

%% filter, epoch and select data:
EEG = pop_eegfilt(EEG, 0, LPF, [], [0]);
EEG = pop_select( EEG,'point',[1 1000] ); % Ausschnitt aus einem EEG erzeugen


% alle EEG latencies einlesen:
a = {EEG.event(1,:).latency}'
cell2mat(a);

%% events hinzufügen:
for i=1:size(EEG.ET.blinks,1)
    if EEG.ET.blinks(i,1)>0 && EEG.ET.blinks(i,2)<EEG.timeAxis(end)
        EEG.event(end+1).type='blink';
        EEG.event(end).latency=EEG.ET.blinks(i,1);
        EEG.event(end).duration=EEG.ET.blinks(i,2)-EEG.ET.blinks(i,1); % das hier nur, wenn allerneuestes EEGLAB 
    end
end
EEG = eeg_checkset(EEG);
% EEG = eeg_checkset(EEG, 'eventconsistency'); % bei älteren Versionen
eeglab redraw;

%% Calculate ERP and Global field power:
% einfaches EKP mit drüber geplotteter GFP:
ERP = mean(EEG.data,3); GFP = std(ERP);
figure; plot(EEG.times,ERP,'k'); 
hold on; plot(EEG.times,GFP,'r','linew',2);

%% calculate min-amplitude, position (index) and latency within time interval
tw1=152; % interval begin (in ms)
tw2=200; % interval end (in ms)
pos1=find(EEG.times==tw1); % find corresponding index 
pos2=find(EEG.times==tw2); % 
min_ch = min(min(erp(:,pos1:pos2))); % find min (respectively max) value
[ch_index timepoint]=find(erp==min_ch); % find corresponding channel and corresponding index value

%%% b) %%%
amp_N170=min_ch;
pos_N170=timepoint; % position (index) of min value 
lat_N170=EEG.times(pos_N170); % latency of min value (in ms);

%% Topoplot für bestimmten Zeitpunkt plotten:
lat = find(EEG.times == 124) % finde Index, wo der Wert 124 steht
figure; topoplot(ERP(:,lat), EEG.chanlocs);

%% How to get selected labels from EEG.event:
a = {EEG.chanlocs([1,2,3]).labels}


%% Handling mehrerer Datensets
% Switche zwischen Datensets hin und her.
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0); 

% lade mehr als ein Datenset rein:
EEG = pop_loadset([filename '_nontarget.set'], cp);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, 0);


%% How to calculate the ica activations:
 tmp = (EEG.icaweights * EEG.icasphere) * EEG.data; % angenommen das EEG.data noch nicht epochiert ist
                                                    % ansonstenm muss man
                                                    % diese Formel für
                                                    % jeden einzelnen
                                                    % Trials ausführen.
                                                    
% the inverse weight matrix EEG.icawinv is equal to
EEG.icawinv = pinv( EEG.icaweights*EEG.icasphere );

% 
EEG.icaweights * EEG.icasphere; % results in an components x channel-weightings matrix.

%% Remove baseline
EEG = pop_rmbase(EEG,[-200 0]);  
