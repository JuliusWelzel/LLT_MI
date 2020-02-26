function [stim] = splitLLTstim(stimuli)
% split LLT stimuli from Braun ea. 2017 for furhter processing
% 4 positions - 1. Laterality [LH, RH, LF, RF] // 2. View [back,palm] // 3. in depth rotation // 4. lateral rotation 
    % stimuli
    for st = 1:length(stimuli)
        stim_pics = stimuli;
        stim_split = strsplit(stim_pics{st},{'_','.'});
        stim(st,:) = stim_split(2:5); % stores 4 properties of stimulus in cell in RT_ALL Position: 1=Limb, 2= orientation, 3=depth rotation, 4 = lateral rotation
    end
    
    stim = stim';
end

