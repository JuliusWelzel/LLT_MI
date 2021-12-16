function [EEG] = LLT_stim( EEG )

% LLT_stim: Function to sort LLT stimuli for certain condition
%     This function helps to differentiate simuli for further analysis.
%     
% INPUT: 
%     OBJ: EEG stricture with EEG.event with LLT stimuli

% OUTPUT:
%     con_ind = condition index see below
%     stim_ind = stimulus index
%     acc = accuracy of LLT response
% 
% 
% Author: (Julius Welzel & Mareike Daeglau, University of Oldenburg, 2018)

%%  find all presentation markers which are NOT pictures and store them

all_trig = {EEG.event(:).type};
not_pic = ismember(all_trig,{'1','2','3','4','empty'}); %all triggers which are not a pic
lat_pic = [EEG.event(~not_pic).latency]; % latency of pictures
type_pic = {EEG.event(~not_pic).type}; %type of stimulus
pic_i = find(~ismember(all_trig,{'1','2','3','4','empty'}));
pic_trig = all_trig(~not_pic);

for ind = 1:length(pic_i)
    EEG.epochs_real(ind).ep_index = ind;
    EEG.epochs_real(ind).trig = EEG.event(pic_i(ind)).type;
    EEG.epochs_real(ind).response = EEG.event(pic_i(ind)+1).type;
end

%find training & experimental pictures
T_pic_all = strfind(all_trig,'T');
E_pic_all = strfind(all_trig,'E');


%% split data per side & extrimity

pic_names = char(type_pic);


%% sort for further analysis

for e = 1:size(pic_trig,2)

% add condition
    if ~isempty(strfind(pic_names(e,:),'lh'));
        EEG.con_ind(e) = 1;
    elseif ~isempty(strfind(pic_names(e,:),'rh'));
        EEG.con_ind(e) = 2;        
    elseif ~isempty(strfind(pic_names(e,:),'lf'));
        EEG.con_ind(e) = 3;;
    else ~isempty(strfind(pic_names(e,:),'rf'));
        EEG.con_ind(e) = 4;
    end

% add accuracy
    if (~isempty(strfind(EEG.epochs_real(e).trig(3),'l')) && strcmp({EEG.epochs_real(e).response},'1'))|(~isempty(strfind(EEG.epochs_real(e).trig(3),'r')) && strcmp({EEG.epochs_real(e).response},'2'));
        EEG.acc (e) = 1;
    else 
        EEG.acc (e) = 0;
    end

% add degrees
% 1 -> medial , 2 -> lateral
    if ~isempty(strfind(EEG.epochs_real(e).trig(end-5),'6'))|~isempty(strfind(EEG.epochs_real(e).trig(end-5),'2'))&&~isempty(strfind(EEG.epochs_real(e).trig(3),'l'))
        EEG.stim_ind(e) = 1;
    elseif ~isempty(strfind(EEG.epochs_real(e).trig(end-5),'4'))|~isempty(strfind(EEG.epochs_real(e).trig(end-5),'0'))&&~isempty(strfind(EEG.epochs_real(e).trig(3),'l'))
        EEG.stim_ind(e) = 2;
    elseif ~isempty(strfind(EEG.epochs_real(e).trig(end-5),'4'))|~isempty(strfind(EEG.epochs_real(e).trig(end-5),'0'))&&~isempty(strfind(EEG.epochs_real(e).trig(3),'r'))
        EEG.stim_ind(e) = 1;
    elseif ~isempty(strfind(EEG.epochs_real(e).trig(end-5),'6'))|~isempty(strfind(EEG.epochs_real(e).trig(end-5),'2'))&&~isempty(strfind(EEG.epochs_real(e).trig(3),'r'))
        EEG.stim_ind(e) = 2;
    end

end

stimuli = pic_trig;
EEG.stim = num2cell(zeros(length(stimuli),4));

for st = 1:length(stimuli)
    stim_split = strsplit(stimuli{st},{'_','.'});
    EEG.stim(st,:) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=lateral rotation, 4 = depth rotation
end


end
