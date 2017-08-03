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

stimConfigs = theData(runNumber).p.stimConfigurations;
finThreshSurr = theData(runNumber).data.finalThresholdSurr;
finThreshNSurr = theData(runNumber).data.finalThresholdNSurr;
[finThreshNSurr finThreshSurr]
cThreshSurr = theData(runNumber).data.cThresholdsSurr;
cThreshNSurr = theData(runNumber).data.cThresholdsNSurr;


figure(1)
plot(cThreshSurr)
title('surround contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])


figure(2)
plot(cThreshNSurr)

title('no surround contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])

figure(3)
bar([finThreshNSurr finThreshSurr])

title('final thresholds')
xlabel('conditions')
ylabel('contrast threshold')
axis square
ylim([0 1])
set(gca, 'XTickLabel', {'no surr' 'surr'})
    

