% Analyze Surround SuppressionV2

clear all
close all

subject = 'Pre-Pilot2_IB';

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

%%
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

% [coll orth ns]
for nRun = 1:length(runNumbers)
    for nConfig = 1:length(configs)
        finThreshQAvgs(1,nConfig,nRun) = mean(allFinThreshQ{nRun}(threshAttIndx(:,nConfig))-fixedStimContrast);
        finThreshQSTD(1,nConfig,nRun) = std(allFinThreshQ{nRun}(threshAttIndx(:,nConfig)));
        finThreshQSTE(1,nConfig,nRun) = finThreshQSTD(1,nConfig,nRun)/nSCAtt;

        finThreshQAvgs(2,nConfig,nRun) = mean(allFinThreshQ{nRun}(threshUnAttIndx(:,nConfig))-fixedStimContrast);
        finThreshQSTD(2,nConfig,nRun) = std(allFinThreshQ{nRun}(threshUnAttIndx(:,nConfig)));
        finThreshQSTE(2,nConfig,nRun) = finThreshQSTD(2,nConfig,nRun)/nSCUnAtt; 
    end
end

%% plots

for nRun = 1:length(runNumbers)
    % final threshold plots
    figure(nRun)
    hold on
    bar(1:3,finThreshQAvgs(:,:,nRun)')
    h = errorbar([1:3;1:3],finThreshQAvgs,finThreshQSTE,'x');
    set(h,'MarkerSize',0.1)
    title(['final thresholds ' subject(end-1:end)])
    xlabel('Condition')
    ylabel('(C_T-C_F)')
    legend('att','unatt')
    axis square
    ylim([0 max(max(finThreshQAvgs(:,:,nRun)))+min(min(finThreshQAvgs(:,:,nRun)))])
    set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
    set(gca, 'XTick', 1:length(configs))
    % staircase plots
    for nStruct = 1:numQStructs
        figure(length(runNumbers)+nRun)
        subplot(4,3,nStruct)
        hold on
        plot(cThreshQ{nRun}(:,nStruct))
        title(qStructNames(nStruct))
        ylabel('C_T')
        xlabel('trial')
        ylim([0 1])
    end
end


