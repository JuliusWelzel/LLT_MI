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
%   'r^2' for biserial correlation coefficients
% chanlocs: channel location information for the topoplot (e.g. typically EEG.chanlocs)
%
% Output:
% sig: a channels x bins matrix yielding the values of the specified
% method.
% fNeu: integer array specifying all frequency bins have been chosen





function [sig,fNeu] = topFreq(data1, data2,f,FOI,method,chanlocs)

    data1 = real(log(data1));
    data2 = real(log(data2));

    % check whether FOI has been specified
    if isempty(FOI)
       FOI(1) = f(1); FOI(2)=f(end);
    end

    %% Select relevant frequency bins
    fNeu = f(f>=FOI(1) & f<=FOI(2));       % only choose interesting frequencies
    idx = find(f>=FOI(1) & f<=FOI(2));              % get corresponding indices

    %% find optimal way to plot all the subplots
    v = ceil(sqrt(length(fNeu)));       % vertical subplots
    h = ceil(sqrt(length(fNeu)));       % horizontal subplots

    %% Topoplot effect strength for selected frequency bins.
    if strcmp(method,'d')
        disp('use strength of effect measure');
        
        % remove all uninteresting frequency bins
        class1 = data1(:,:,idx);         % channels x trials x selected frequency bins
        class2 = data2(:,:,idx);         % channels x trials x selected frequencys bins

        % calculate strength of difference
        mean1 = squeeze(mean(class1,2)); % channels x frequency bins
        mean2 = squeeze(mean(class2,2)); % channels x freqzency bins
        std1 = squeeze(std(class1,0,2)); % channels x frequency bins
        std2 = squeeze(std(class2,0,2)); % channels x frequency bins
        effect = (mean1 - mean2)./sqrt((std1.^2+std2.^2)./2); % channels x frequency bins

        % Frontal electrodes
        effect([53,45,8,52,16,44,29,30,9,46],:) = 0;
        
        % Posterior electrodes
        effect([49,57,56],:) = 0;
        
        % calculate overall power
        oP = mean(abs(effect),1); % 1 x bins
        
        % find appropriate common scaling for all topoplots
        tmpMin = min(effect(:)); tmpMax = max(effect(:));
        %tmpMin = mean(effect(:))- 1.5*std(effect(:)); tmpMax = mean(effect(:))+1.5*std(effect(:))%max(effect(:));

        % plot amplitude differences 
        %figure(21);
        for i=1:length(fNeu)

            % update topoplot
            subplot(v+1,h,i)

            if i == length(fNeu)
                topoplot(effect(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
                hc = colorbar;
               % colorbar([0.1 0.1  0.7  0.8]);
                set(hc,'Position',[0.8 0.1 0.05 0.1])
            else
                topoplot(effect(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
            end

            % update corresponding frequency
            title([num2str(fNeu(i)) ' Hz']);    
        end

        subplot(v+1,h,v*h+1:(v+1)*h-1);
        plot(fNeu,oP);
        title('Most informative frequency bins','fontsize',15); xlabel('Frequency (Hz)'); ylabel('|d|');
   
        
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.45,0.97,'Effect strength (cohens d): class1>class2','fontsize',15)
        sig = effect;
                
        
    end

    %%% Differences as visualised in the significance level 
    if strcmp(method,'p')
        disp('use level of significance measure');
        
        % select and reduce data to intersting frequency bins
        sig1 = data1(:,:,idx);
        sig2 = data2(:,:,idx);

        % bring data into write format for ttest()
        sig1 = permute(sig1,[2 1 3]);
        sig2 = permute(sig2,[2 1 3]);

        % calc ttest
        [tmp,p] = ttest(sig1,sig2); p = squeeze(p); % p --> channels x frequency bins

        % find appropriate common scaling for all topoplots
        tmpMin = 0; tmpMax = 0.2;

        %figure(21);
        for i=1:length(fNeu)

            % update topoplot
            subplot(v,h,i)


            if i == length(fNeu)
                topoplot(p(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
                hc = colorbar;
               % colorbar([0.1 0.1  0.7  0.8]);
                set(hc,'Position',[0.8 0.1 0.05 0.1])
            else
                topoplot(p(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
            end

            % update corresponding frequency
            title([num2str(fNeu(i)) ' Hz']);

        end
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.45,0.97,'Level of significance of amplitude differences','fontsize',15)
        sig = p;

    end


    %% Topoplot amplitude differences for selected frequency bins.
    if strcmp(method,'diff')
        disp('use amplitude difference measure');
        
        % average across trials and squeeze matrix
        class1 = squeeze(mean(data1,2));    %  channels x frequency bins
        class2 = squeeze(mean(data2,2));     %  channels x frequency bins

        % remove all uninteresting frequency bins
        class1 = class1(:,idx);
        class2 = class2(:,idx);

        % calculate amplitude differences between conditions:
        ampDiff = class1 - class2;            % channels x frequency bins


        % logarithmize for better visualisation:
        %class1 = log(class1); class2 = log(class2);
        %ampDiff = log10(ampDiff);

        % find appropriate common scaling for all topoplots
        tmpMin = min(ampDiff(:)); tmpMax = max(ampDiff(:));

        % plot amplitude differences 
        %figure(21);
        for i=1:length(fNeu)

            subplot(v,h,i) 
            if i == length(fNeu)
                topoplot(ampDiff(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
                h = colorbar;
                set(h,'Position',[0.8 0.1 0.05 0.1])
            else
                topoplot(ampDiff(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
            end

            % update corresponding frequency
            title([num2str(fNeu(i)) ' Hz']);

        end
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.45,0.97,'Amplitude differences: class1-class2','fontsize',15)
        sig = ampDiff;
    end
    
    %% Topoplot point-biseral correlations
    if strcmp(method,'r')
        disp('use Point-biserial correlation coefficient');

        % select and reduce data to interesting frequency bins
        sig1 = data1(:,:,idx);
        sig2 = data2(:,:,idx);
        
        c0 = ones(size(sig1));
        c1 = zeros(size(sig2));
        
        myCorr = [];
        for c=1:size(sig1,1)        % iterate through channels
            for b=1:size(sig1,3)    % iterate through bins
                tmp = corrcoef([sig1(c,:,b), sig2(c,:,b)],[c0(c,:,b), c1(c,:,b)]);
  
                myCorr(c,b) = tmp(2); 
            end
        end
        
        % square to r^2
        myCorr = myCorr.^2;       
                
        % calculate overall power
        %OP = mean(abs(myCorr),1); % 1 x bins
        
        % calulcate power for electrode sites C3,Cz,C4
        chanLab = {chanlocs.labels};
        C3 = find(strcmp(chanLab,'E55')); C3 = myCorr(C3,:); 
        CZ = find(strcmp(chanLab,'E01')); CZ = myCorr(CZ,:);
        C4 = find(strcmp(chanLab,'E47')); C4 = myCorr(C4,:);
              
        
        % find appropriate common scaling for all topoplots
        tmpMin = min(myCorr(:)); tmpMax = max(myCorr(:));

        
        % plot amplitude differences 
        %figure(21);
        for i=1:length(fNeu)

            subplot(v+1,h,i) 
            if i == length(fNeu)
                topoplot(myCorr(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
                hc = colorbar;
                set(hc,'Position',[0.8 0.1 0.05 0.1])
            else
                topoplot(myCorr(:,i),chanlocs, 'electrodes','off','maplimits',[tmpMin tmpMax]);
            end

            % update corresponding frequency
            title([num2str(fNeu(i)) ' Hz']);

        end
        
        subplot(v+1,h,v*h+1:(v+1)*h-1);
        plot(fNeu,[C3; CZ; C4]);
        legend('C3', 'Cz', 'C4');
        title('Most informative frequency bins','fontsize',15); xlabel('Frequency (Hz)'); ylabel('r^2');
        
        axes('Position',[0 0 1 1],'Visible','off');
        text(0.45,0.97,'Point-biserial correlation (r^2)','fontsize',15)
        sig = myCorr;
    end
    
    
    