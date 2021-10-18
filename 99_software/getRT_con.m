function [ sixty onetwnty twoforty threehun ] = getRT_con( RT_ALL, deg)
%get_RT_con Function to get RT of different stimuli for different
%conditions

% Extracting RT of RT_ALL matrix of nic_LLT dataset
% Depenging on the question, subsets of data are
% required. This function helps to collect important data suiting your
% research question

% INPUT:
% RT_ALL = Matrix of all stimuli of all participants
% deg = Rotation degree in depth of stimuli
% list_sub = number of subjects to average over


% Author: Julius Welzel, Uni Oldenburg, 2018

%% Extract degree degree roation horizontally
degree = string(deg);
c = 1;

for sub = 1:length(RT_ALL)
    for st = 10:length(RT_ALL(sub).stim)
        if sum(strcmp(RT_ALL(sub).stim{st,3},degree))>=1  % only take degrees horizontally which are wanted
            RT_DS(c).deg = str2num(RT_ALL(sub).stim{st,4}); %get the degree of rotation
            RT_DS(c).ID = str2num(RT_ALL(sub).ID(5:6)); %get the Subject ID
            RT_DS(c).limb = RT_ALL(sub).stim{st,1}; %get limb
            RT_DS(c).mod = RT_ALL(sub).con_ind(st); %get the degree of lateralt = 1 /medial = 2
            RT_DS(c).RT = RT_ALL(sub).SO_ms(st); %get RT 
            RT_DS(c).orient = RT_ALL(sub).stim{st,2}; % get orientation of stimulus back/palm
            RT_DS(c).acc = RT_ALL(sub).acc(st);

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
        
sixty.lh_o = [sixty_lh_o{:}];
% sixty_lh_o(isoutlier(sixty_lh_o)) = [];
sixty.lh_y = [sixty_lh_y{:}];
sixty.rh_o = [sixty_rh_o{:}];
sixty.rh_y = [sixty_rh_y{:}];

onetwnty.lh_o = [onet_lh_o{:}];
% onet_lh_o(isoutlier(onet_lh_o)) = [];
onetwnty.lh_y = [onet_lh_y{:}];
onetwnty.rh_o = [onet_rh_o{:}];
onetwnty.rh_y = [onet_rh_y{:}];

twoforty.lh_o = [twofor_lh_o{:}];
% twofor_lh_o(isoutlier(twofor_lh_o)) = [];
twoforty.lh_y = [twofor_lh_y{:}];
twoforty.rh_o = [twofor_rh_o{:}];
twoforty.rh_y = [twofor_rh_y{:}];

threehun.lh_o = [three_lh_o{:}];
% three_lh_o(isoutlier(three_lh_o)) = [];
threehun.lh_y = [three_lh_y{:}];
threehun.rh_o = [three_rh_o{:}];
threehun.rh_y = [three_rh_y{:}];

end

