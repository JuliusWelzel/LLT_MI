function [extremity condition angle]  = LLTepochinfo(ep_stim)

% extrimity
    if contains(ep_stim(1),'h')
        extremity = "hand";
    elseif contains(ep_stim(1),'f')
        extremity = "foot";
    end
    
% condition
    if contains(ep_stim(1),'l') & contains(ep_stim(4),{'240','300'})|...
        contains(ep_stim(1),'r') & contains(ep_stim(4),{'60','120'});
        condition = "lateral";
    elseif contains(ep_stim(1),'l') & contains(ep_stim(4),{'60','120'})|...
        contains(ep_stim(1),'r') & contains(ep_stim(4),{'240','300'});
        condition = "medial";
    end
    
% angle
    if strcmp(ep_stim(4),'60')| strcmp(ep_stim(4),'300')
        angle = "60";
    elseif strcmp(ep_stim(4),'120')| strcmp(ep_stim(4),'240')
        angle = "120";
    end

end

