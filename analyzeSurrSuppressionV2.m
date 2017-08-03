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
finThreshSurrCu = theData(runNumber).data.finalThresholdSurrCu;
finThreshNSurrCu = theData(runNumber).data.finalThresholdNSurrCu;
finThreshSurr = theData(runNumber).data.finalThresholdSurr;
finThreshNSurr = theData(runNumber).data.finalThresholdNSurr;

cThreshSurrCu = theData(runNumber).data.cThresholdsSurrCu;
cThreshNSurrCu = theData(runNumber).data.cThresholdsNSurrCu;
cThreshSurr = theData(runNumber).data.cThresholdsSurr;
cThreshNSurr = theData(runNumber).data.cThresholdsNSurr;

%% plots

figure(1)
plot(cThreshSurrCu)
title('surround with cue contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])


figure(2)
plot(cThreshNSurrCu)
title('no surround with cue contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])

figure(3)
plot(cThreshSurr)
title('surround no cue contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])


figure(4)
plot(cThreshNSurr)
title('no surround no cue contrast thresholds')
xlabel('trial')
ylabel('contrast threshold')
axis square
ylim([0 1])

figure(5)
bar([finThreshNSurrCu finThreshSurrCu finThreshNSurr finThreshSurr])
title('final thresholds')
xlabel('conditions')
ylabel('contrast threshold')
axis square
ylim([0 1])
set(gca, 'XTickLabel', {'no surr cued' 'surr cued' 'no surr co cue' 'surr no cue'})
    

