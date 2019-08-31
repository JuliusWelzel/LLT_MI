% Topoplots the differences between two epoched datasets for specific
% frequency bins
% 
% Input:
% data1: 3D array of power values for 1st class (channels x trials x frequency bins)
% data2: 3D array of power values for 2nd class (channels x trials x frequency bins)
% f: integer array specifying all frequency bins (must conver to the number of frequency bins in data1,2)
% FOI: either [] or 2-elements array specifying the subset of Frequencies
% of interest
% method: a string specifying the method: 
%   'p' for level of significance,
%   'd' for strength of effect, 
%   'diff' for amplitude differences and 
%   'corr' for biserial correlation coefficients
% chanlocs: channel location information for the topoplot (e.g. typically EEG.chanlocs)
%
% Output:
% sig: a channels x bins matrix yielding the values of the specified
% method.
%




function sig = topFreq2(data,f,chanlocs)

    x = mean(data(:))+std(data(:));
    data(data<x) = 0; 
    % plot amplitude differences 
    %figure(21);
    for i=1:length(fNeu)

        subplot(v,h,i) 
        if i == length(fNeu)
            topoplot(data(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
            h = colorbar;
            set(h,'Position',[0.8 0.1 0.05 0.1])
        else
            topoplot(data(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
        end

        % update corresponding frequency
        title([num2str(f(i)) ' Hz']);

    end
    axes('Position',[0 0 1 1],'Visible','off');
    text(0.45,0.97,'Signal stength','fontsize',15)
    sig = myCorr;

    
    