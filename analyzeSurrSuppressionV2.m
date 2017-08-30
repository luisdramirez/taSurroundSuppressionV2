% Analyze Surround SuppressionV2

clear all
close all

subject = 'Pilot_LV';

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppressionV2_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppressionV2_', subject, '.mat']);
    runNumbers = 1:length(theData);
else
    error('Data file does not exist.')
end

%% prepare data 
for nRun = 1:length(runNumbers)
    data{nRun} = theData(nRun).data;
    p{nRun} = theData(nRun).p;
end

stimConfigs = p{1}.stimConfigurations;
fixedStimContrast = p{1}.fixedStimContrast;
numQStructs = p{1}.numQStructures;
qStructNames = {'collCued1' 'orthCued1' 'nsCued1' 'collCued2' 'orthCued2' 'nsCued2' 'collCued3' 'orthCued3' 'nsCued3' 'collUnCued1' 'orthUnCued1' 'nsUnCued1'}; %theData(runNumber).p.qStructuresNames;

for nRun = 1:length(runNumbers)
    allFinThreshQ{nRun} = data{nRun}.finalThresholdQ;
    cThreshQ{nRun} = data{nRun}.cThresholdsQ;
end

nSCAtt = 3;
nSCUnAtt = 1;
configs = {'coll' 'orth' 'ns'};
cues = {'Cued' 'UnCued'};

threshAttIndx = nan(nSCAtt,length(configs));
threshUnAttIndx = nan(nSCUnAtt,length(configs));

% determine which final thresholds to take as averages for each
% combination
for nConfig = 1:length(configs)
    for nSC = 1:nSCAtt
        currConfig = [configs{nConfig} cues{1} num2str(nSC)];
        for nQStruct = 1:length(qStructNames)
            currStruct = qStructNames{nQStruct};
            if strcmp(currConfig,currStruct)
                threshAttIndx(nSC,nConfig) = nQStruct;
            end
        end
    end
    for nSC = 1:nSCAtt
        currConfig = [configs{nConfig} cues{2} num2str(nSC)];
        for nQStruct = 1:length(qStructNames)
            currStruct = qStructNames{nQStruct};
            if strcmp(currConfig,currStruct)
                threshUnAttIndx(nSC,nConfig) = nQStruct;
            end
        end
    end
end

% final threshold averages, std, and ste [coll orth ns]
for nConfig = 1:length(configs)
    finThreshQAvgs(1,nConfig) = mean(allFinThreshQ(threshAttIndx(:,nConfig))-fixedStimContrast);
    finThreshQSTD(1,nConfig) = std(allFinThreshQ(threshAttIndx(:,nConfig)));
    finThreshQSTE(1,nConfig) = finThreshQSTD(1,nConfig)/nSCAtt;
    
    finThreshQAvgs(2,nConfig) = mean(allFinThreshQ(threshUnAttIndx(:,nConfig))-fixedStimContrast);
    finThreshQSTD(2,nConfig) = std(allFinThreshQ(threshUnAttIndx(:,nConfig)));
    finThreshQSTE(2,nConfig) = finThreshQSTD(2,nConfig)/nSCUnAtt; 
end

% determine average thresholds per condition, per run [(att unatt) x (coll orth ns)]
for nRun = 1:length(runNumbers)
    for nConfig = 1:length(configs)
        finThreshQAvgs(1,nConfig,nRun) = mean(allFinThreshQ{nRun}(threshAttIndx(:,nConfig)))-fixedStimContrast;
        finThreshQSTD(1,nConfig,nRun) = std(allFinThreshQ{nRun}(threshAttIndx(:,nConfig)));
        finThreshQSTE(1,nConfig,nRun) = finThreshQSTD(1,nConfig,nRun)/nSCAtt;

        finThreshQAvgs(2,nConfig,nRun) = mean(allFinThreshQ{nRun}(threshUnAttIndx(:,nConfig)))-fixedStimContrast;
        finThreshQSTD(2,nConfig,nRun) = std(allFinThreshQ{nRun}(threshUnAttIndx(:,nConfig)));
        finThreshQSTE(2,nConfig,nRun) = finThreshQSTD(2,nConfig,nRun)/nSCUnAtt; 
    end
end

% determine percent difference between att and unatt per condition [run x (coll orth ns)]
for nRun = 1:length(runNumbers)
    for nConfig = 1:length(configs)
        percDiffFinThreshQAvgs(nRun,nConfig) = (finThreshQAvgs(2,nConfig,nRun)-finThreshQAvgs(1,nConfig,nRun))/finThreshQAvgs(1,nConfig,nRun);
    end
end

% determine average thresholds per condition across runs [(att unatt) x (coll orth ns)]
if length(runNumbers) > 1
    allFinThreshQAvgs = mean(finThreshQAvgs,3);
    allFinThreshQSTD = std(finThreshQAvgs,0,3);
    allFinThreshQSTE = allFinThreshQSTD/length(runNumbers);
    for nConfig = 1:length(configs)
        allPercDiffFinThreshQAvgs(nConfig) = (allFinThreshQAvgs(2,nConfig)-allFinThreshQAvgs(1,nConfig))/allFinThreshQAvgs(1,nConfig);
    end
end

%% plots

% Plot each run
for nRun = 1:length(runNumbers)
    % final threshold plots (difference)
    figure
    bar(1:3,finThreshQAvgs(:,:,nRun)')
    hold on
    h = errorbar([1:3;1:3],finThreshQAvgs(:,:,nRun),finThreshQSTE(:,:,nRun),'x');
    set(h,'MarkerSize',0.1)
    title(['run: ' num2str(nRun) ' final thresholds ' subject(end-1:end)])
    xlabel('Condition')
    ylabel('(C_T-C_F)')
    legend('att','unatt')
    axis square
    ylim([0 max(max(finThreshQAvgs(:,:,nRun)))+min(min(finThreshQAvgs(:,:,nRun)))])
    set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
    set(gca, 'XTick', 1:length(configs))
    hold off
    % staircase plots
    figure
    for nStruct = 1:numQStructs
        subplot(4,3,nStruct)
        plot(cThreshQ{nRun}(:,nStruct))        
        title(qStructNames(nStruct))
        ylabel('C_T')
        xlabel('trial')
        ylim([0 max(cThreshQ{nRun}(:,nStruct))+min(cThreshQ{nRun}(:,nStruct))])
    end
    % final threshold plots (percent difference)
    figure
    hold on
    bar(1:3,percDiffFinThreshQAvgs(nRun,:)')
    title(['run: ' num2str(nRun) ' perc diff final thresholds ' subject(end-1:end)])
    xlabel('Condition')
    ylabel('% diff')
    axis square
    ylim([0 max(max(percDiffFinThreshQAvgs(nRun,:)))+min(min(percDiffFinThreshQAvgs(nRun,:)))])
    set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
    set(gca, 'XTick', 1:length(configs))   
    hold off
end

% Plot all runs
if length(runNumbers) > 1
    % final threshold plots (difference)
    figure
    hold on
    bar(1:3,allFinThreshQAvgs')
    h2 = errorbar([1:3;1:3],allFinThreshQAvgs,allFinThreshQSTE,'x');
    set(h2,'MarkerSize',0.1)
    title(['all runs (' num2str(runNumbers(end)) ') final thresholds ' subject(end-1:end)])
    xlabel('Condition')
    ylabel('(C_T-C_F)')
    legend('att','unatt')
    axis square
    ylim([0 max(max(allFinThreshQAvgs))+min(min(allFinThreshQAvgs))])
    set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
    set(gca, 'XTick', 1:length(configs))
    
    % final threshold plots (percent difference)
    figure
    hold on
    bar(1:3,percDiffFinThreshQAvgs(nRun,:)')
    title(['all runs (' num2str(runNumbers(end)) ') perc diff final thresholds ' subject(end-1:end)])
    xlabel('Condition')
    ylabel('% diff')
    axis square
    ylim([0 max(max(percDiffFinThreshQAvgs(nRun,:)))+min(min(percDiffFinThreshQAvgs(nRun,:)))])
    set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
    set(gca, 'XTick', 1:length(configs))   
    hold off
end

