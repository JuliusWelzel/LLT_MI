% Plots the spectral differences between two classes on a specific ROI
% Input:
% data1: 3D array of power values for 1st class (channels x trials x frequency bins)
% data2: 3D array of power values for 2nd class (channels x trials x frequency bins)
% f: integer array specifying all frequency bins (must conver to the number of frequency bins in data1,2)
% FOI: either [] or 2-elements array specifying the subset of Frequencies
% of interest
% chanlocs: channel location information for the topoplot (e.g. typically EEG.chanlocs)
% method: string 'd' for stength of effect or 'corr' for biserial
% correlation coefficients

function plotROI(data1,data2,f,FOI,ROI,chanlocs,method)

    % check whether FOI has been specified
    if isempty(FOI)
       FOI(1) = f(1); FOI(2)=f(end);
    end

    % select FOI
    fNeu = f(f>=FOI(1) & f<=FOI(2));     % only choose interesting frequencies
    idx = find(f>=FOI(1) & f<=FOI(2));     % get corresponding indices
    
    if strcmp(method,'d')
    
        % select relevant ROI'S
        class1 = squeeze(mean(data1(ROI,:,idx),1));  % calculate amplitude values for ROI
                                                     % trials x frequency bins

        class2 = squeeze(mean(data2(ROI,:,idx),1));  % calculate amplitude values for ROI
                                                          % trials x frequency bins

        % Calculate strength of effect:
        effect = (mean(class1) - mean(class2))./sqrt((std(class1).^2+std(class2).^2)./2);

        % Logarithmize for better visualisation                                              
        class1 = log10(class1); class2 = log10(class2);    % trials x frequency bins

        % calculate STD across trials
        class1Std = std(class1);                           % 1 x frequency bins
        class2Std = std(class1);                           % 1 x frequency bins

        % Prepare topoplot for visualizing ROI
        topoROI = zeros(1,size(data1,1)); topoROI(ROI) = 1;

        %figure(22);
        subplot(4,1,[1 2 3]);
        shadedErrorBar(fNeu,mean(class1,1),class1Std, 'r',1); hold on;
        shadedErrorBar(fNeu,mean(class2,1),class2Std, 'b',1); hold on;
        plot(fNeu,mean(class1,1),'r','linewidth',3);
        plot(fNeu,mean(class2,1),'b','linewidth',3);
        ylabel('Amplitude');
        %xlabel('Frequency');
        title(['Power spectral density (channels ' num2str(ROI) ')'],'fontsize',15);

        subplot(4,1,4);
        plot(fNeu,effect);
        hline(0,'r');
        ylabel('Cohens d');
        xlabel('Frequency (Hz)');
        title('Strength of effect: class1 (r) >class2 (b)','fontsize',15);

        if ~isempty(chanlocs)
            axes('position',[.6  .6  .4  .4]);
            topoplot(topoROI,chanlocs,'electrodes','numbers');
        end
    % Not yet ready to be used ( use cohens d instead)!!!
    elseif strcmp(method,'r')
         
        % select relevant ROI'S
        class1 = squeeze(mean(data1(ROI,:,idx),1));  % calculate amplitude values for ROI
                                                     % trials x frequency bins

        class2 = squeeze(mean(data2(ROI,:,idx),1));  % calculate amplitude values for ROI
                                                     % trials x frequency bins
        
        % Calculate biserial correlation
        corr = [];
        c1 = ones(1,size(class1,1));
        c2 = zeros(1,size(class2,1));
        
        keyboard;
        for b=1:size(class1,2)            
            tmpCorr = corrcoef([class1(:,b);class1(:,b)],[c1, c2]');
            corr(b) = tmpCorr(2);
        end
        
        % Logarithmize for better visualisation                                              
        class1 = log10(class1); class2 = log10(class2);    % trials x frequency bins

        % calculate STD across trials
        class1Std = std(class1);                           % 1 x frequency bins
        class2Std = std(class1);                           % 1 x frequency bins

        % Prepare topoplot for visualizing ROI
        topoROI = zeros(1,size(data1,1)); topoROI(ROI) = 1;
        
        % Visualise
        subplot(4,1,[1 2 3]);
        shadedErrorBar(fNeu,mean(class1,1),class1Std, 'r',1); hold on;
        shadedErrorBar(fNeu,mean(class2,1),class2Std, 'b',1); hold on;
        plot(fNeu,mean(class1,1),'r','linewidth',3);
        plot(fNeu,mean(class2,1),'b','linewidth',3);
        ylabel('Amplitude');
        %xlabel('Frequency');
        title(['Power spectral density (channels ' num2str(ROI) ')'],'fontsize',15);

        subplot(4,1,4);
        plot(fNeu,corr);
        hline(0,'r');
        ylabel('r');
        xlabel('Frequency (Hz)');
        title('Biserial correlation','fontsize',15);

        if ~isempty(chanlocs)
            axes('position',[.6  .6  .4  .4]);
            topoplot(topoROI,chanlocs,'electrodes','numbers');
        end
        
        
        
        
        
        
        
    end        
    