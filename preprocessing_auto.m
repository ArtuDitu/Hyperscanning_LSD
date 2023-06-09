%% Load data
cd /Users/artur/Dropbox/Projects/Hyperscanning_Maastricht
addpath(genpath('D:\Dropbox\Projects\Hyperscanning_Maastricht'))
eeglabpath = fileparts(which('eeglab.m'));
% load eeglab
eeglab;

% load data
EEG = pop_loadbv('/Users/artur/Dropbox/Projects/Hyperscanning_Maastricht/Artur 2/', 'C1_A1_AFF.vhdr')
    
    
    
    
    eeg_sub1 = pop_select(EEG, 'channel', {EEG.chanlocs(1:35).labels});
    eeg_sub2 = pop_select(EEG, 'channel', {EEG.chanlocs(36:70).labels});
    eeg_sub2.chanlocs = eeg_sub1.chanlocs;    
    
    % add events to both datasets (300 events separated by 1 second
    eeg_sub1_chanlocs = eeg_sub1.chanlocs;
    eeg_sub1.data(36,[1000:500:500*floor(size(eeg_sub1.data,2)/500)]) = 1; % simulating a stimulus onset every second
    eeg_sub1 = pop_chanevent(eeg_sub1, 36,'edge','leading','edgelen',0,'duration','on','nbtype',1);
    eeg_sub1.chanlocs = eeg_sub1_chanlocs;
    
    eeg_sub2_chanlocs = eeg_sub2.chanlocs;
    eeg_sub2.data(36,[1000:500:500*floor(size(eeg_sub1.data,2)/500)]) = 1; % simulating a stimulus onset every second
    eeg_sub2 = pop_chanevent(eeg_sub2, 36,'edge','leading','edgelen',0,'duration','on','nbtype',1);
    eeg_sub2.chanlocs = eeg_sub2_chanlocs;
    
    % parameters
    low_pass = 100;
    high_pass = .1;
    srate = 500;
    power_line = 50;
    
    for sub = 1:2
        
        if sub == 1
            EEG = eeg_sub1;
        end
        if sub == 2
            EEG = eeg_sub2;
        end

        EEG.setname

    EEG = pop_chanedit(EEG, 'lookup',fullfile(eeglabpath,'plugins/dipfit/standard_BESA/standard-10-5-cap385.elp'));
        
        % HIGH- AND LOW-PASS FILTERING
        EEG = pop_eegfiltnew(EEG, high_pass, []); % 0.1 is the lower edge
        EEG = pop_eegfiltnew(EEG, [], low_pass); % 100 is the upper edge
        % remove line noise with zapline
        EEG = clean_data_with_zapline_plus_eeglab_wrapper(EEG);
        

    full_chanlocs = EEG.chanlocs;


        
        count = 1;
        
        for event_id = 1:numel(EEG.event)
            if isequal(EEG.event(event_id).type, 'chan36')
                EEG.event(event_id).epoch_id = count;
                count = count + 1;
            else
                EEG.event(event_id).epoch_id = 0;
            end
            
        end
        
        
        % ASR
        EEG = clean_artifacts(EEG);
        
        % save removed channels
        removed_channels = ~ismember({full_chanlocs.labels},{EEG.chanlocs.labels});
        EEG.removed_channels = {full_chanlocs(removed_channels).labels};
        
        % high pass 2 Hz for data used for ICA calculations
        eeg_tmp = pop_eegfiltnew(EEG, 2, []);   % highpass  2 Hz to not include slow drifts
        % create amica folder
        cd D:\Dropbox\Projects\Hyperscanning_Maastricht\amicaFolder
        mkdir(sprintf('ICA_amica_%s_%d',EEG.setname, sub))
        outDir = what(sprintf('ICA_amica_%s_%d',EEG.setname,sub));
        %run ICA
        dataRank = rank(double(eeg_tmp.data'));
        runamica15(eeg_tmp.data, 'num_chans', eeg_tmp.nbchan,...
            'outdir', outDir.path,...
            'pcakeep', dataRank, 'num_models', 1,...
            'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
        
        %load ICA results
        mod = loadmodout15(outDir.path);
        
        %apply ICA weights to data
        EEG.icasphere = mod.S;
        EEG.icaweights = mod.W;
        EEG = eeg_checkset(EEG);
        
        % calculate iclabel classification
        EEG = iclabel(EEG);
        
        % list components that should be rejected
        components_to_remove = [];
        number_components = size(EEG.icaact,1);
        
        for component = 1:number_components
            % muscle
            if EEG.etc.ic_classification.ICLabel.classifications(component,2) > .80
                components_to_remove = [components_to_remove component];
            end
            % eye
            if EEG.etc.ic_classification.ICLabel.classifications(component,3) > .9
                components_to_remove = [components_to_remove component];
            end
            % heart
            if EEG.etc.ic_classification.ICLabel.classifications(component,4) > .9
                components_to_remove = [components_to_remove component];
            end
            % line noise
            if EEG.etc.ic_classification.ICLabel.classifications(component,5) > .9
                components_to_remove = [components_to_remove component];
            end
            % channel noise
            if EEG.etc.ic_classification.ICLabel.classifications(component,6) > .9
                components_to_remove = [components_to_remove component];
            end
        end
        
        % remove components
        EEG = pop_subcomp( EEG, components_to_remove, 0);
        % save removed components in struct
        EEG.removed_components = components_to_remove;
        
        % interpolate removed channels
        EEG = pop_interp(EEG, full_chanlocs,'spherical');
        
        EEG = pop_editset(EEG, 'setname', sprintf('cleaned_%s_sub_%d', EEG.setname, sub));
        EEG = pop_saveset(EEG, 'filename', EEG.setname);
        if sub == 1
            eeg_sub1 = EEG;
        elseif sub == 2
            eeg_sub2 = EEG;
        end
    end
    
    %EPOCH %%%
    window = [0 .999];
    eeg_sub1 = pop_epoch( eeg_sub1, {'chan25'}, window, 'epochinfo', 'yes');
    eeg_sub2 = pop_epoch( eeg_sub2, {'chan25'}, window, 'epochinfo', 'yes');
    
    % finds trials that are present in both subjects
    [val1 pos1] = intersect([eeg_sub1.event.epoch_id], [eeg_sub2.event.epoch_id]);
    [val2 pos2] = intersect([eeg_sub2.event.epoch_id], [eeg_sub1.event.epoch_id]);
    
    % remove trials that are not in both datasets
    eeg_sub1.data = eeg_sub1.data(:,:,pos1);
    eeg_sub1.event = eeg_sub1.event(pos1);
    eeg_sub1.epoch = eeg_sub1.epoch(pos1);
    
    eeg_sub2.data = eeg_sub2.data(:,:,pos2);
    eeg_sub2.event = eeg_sub2.event(pos2);
    eeg_sub2.epoch = eeg_sub2.epoch(pos2);
    
    % create new latencies for concanated data
    latency_sub1 = [];
    count = -499;
    for i = 1:size(eeg_sub1.event,2)
        count = count + 500;
        latency_sub1(i) = count;
    end
    latency_sub2 = [];
    count = -499;
    for i = 1:size(eeg_sub2.event,2)
        count = count + 500;
        latency_sub2(i) = count;
    end
    
    latency_sub1 = num2cell(latency_sub1');
    latency_sub2 = num2cell(latency_sub2');
    
    % concanate data
    eeg_sub1.data = reshape(eeg_sub1.data, 24, []);
    eeg_sub2.data = reshape(eeg_sub2.data, 24, []);
    
    % add new latencies
    [eeg_sub1.event.latency] = latency_sub1{:};
    [eeg_sub2.event.latency] = latency_sub2{:};
    
    eeg_sub1.trials = size(eeg_sub1.data,2)/500;
    eeg_sub2.trials = size(eeg_sub2.data,2)/500;
    
    
    % i commented out line 495 in the script 'po_editset': '%EEGOUT =
    % eeg_checkset(EEGOUT);' because checkset was not allowing for manual
    % concatation of data - change that if you run the script but remember to
    % change it back afterwards.
    
    % and line 132 in pop_saveset
    
    eeg_sub1 = pop_editset(eeg_sub1, 'setname', sprintf('hyper_%s', eeg_sub1.setname));
    eeg_sub1 = pop_saveset(eeg_sub1, 'filename', eeg_sub1.setname);
    eeg_sub2 = pop_editset(eeg_sub2, 'setname', sprintf('hyper_%s', eeg_sub2.setname));
    eeg_sub2 = pop_saveset(eeg_sub2, 'filename', eeg_sub2.setname);

