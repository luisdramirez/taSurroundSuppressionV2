% Analyze Surround SuppressionV2


subject = 'Pre-Pilot_LR';

plotData = 'Yes';

expDir = pwd;
dataDir = 'data';
cd(dataDir)

if exist(['vTA_surrSuppression_', subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppression_', subject, '.mat']);
    runNumber = length(theData);
else
    error('Data file does not exist.')
end

%%
stimConfigs = theData(runNumber).p.stimConfigurations;
allFinThreshQ = theData(runNumber).data.finalThresholdQ;
cThreshQ = theData(runNumber).data.cThresholdsQ;
numQStructs = theData(runNumber).p.numQStructures;
qStructNames = {'collCued1' 'orthCued1' 'nsCued1' 'collCued2' 'orthCued2' 'nsCued2' 'collCued3' 'orthCued3' 'nsCued3' 'collUnCued1' 'orthUnCued1' 'nsUnCued1'}; %theData(runNumber).p.qStructuresNames;

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
for nConfig = 1:length(configs)
    finThreshQAvgs(1,nConfig) = mean(allFinThreshQ(threshAttIndx(:,nConfig)));
    finThreshQSTD(1,nConfig) = std(allFinThreshQ(threshAttIndx(:,nConfig)));
    finThreshQSTE(1,nConfig) = finThreshQSTD(1,nConfig)/nSCAtt;
    
    finThreshQAvgs(2,nConfig) = mean(allFinThreshQ(threshUnAttIndx(:,nConfig)));
%     finThreshQSTD(2,nConfig) = std(allFinThreshQ(threshUnAttIndx(:,nConfig)));
%     finThreshQSTE(2,nConfig) = finThreshQSTD(2,nConfig)/nSCUnAtt; 
end


%% plots

figure(1)
hold on
bar(finThreshQAvgs')
% errorbar(finThreshQAvgs(1,:)',finThreshQSTE(1,:)','.')
title('final thresholds')
xlabel('Condition')
ylabel('C_T')
legend('att','unatt')
axis square
ylim([0 1])
set(gca, 'XTickLabel', {'coll' 'orth' 'ns'})
set(gca, 'XTick', [1:length(configs)])

figure(2)
for nStruct = 1:numQStructs
subplot(4,3,nStruct)
plot(cThreshQ(:,nStruct))
title(qStructNames(nStruct))
ylabel('C_T')
xlabel('trial')
ylim([0 1])
end
