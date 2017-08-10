
% How does temporal attention influence tuned normalization?

%% PREPARE AND COLLECT INFO

commandwindow
HideCursor;
echo off
clear all
close all
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 0);

% Subject name
p.subject = 'Pre-Pilot_LR';
% Trial Events Parameters
p.cueValidity = 0.75;
[p.numAttTrialsPerComb, p.minNumBlocks] = rat(p.cueValidity);
p.repetitions = 20; %20 for at least 40 trials per staircase
p.numBlocks = p.minNumBlocks*p.repetitions; 
p.numQStructures = 12;

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames, ~] = GetKeyboardIndices;
% deviceString = 'Corsair Corsair K95W Gaming Keyboard';
deviceString = 'Apple Inc. Apple Keyboard';
% deviceString = 'Apple Keyboard';
% deviceString = 'CHICONY USB Keyboard';
% deviceString = 'Apple Internal Keyboard / Trackpad';

for i = 1:length(productNames)
    if strcmp(productNames{i}, deviceString)
        deviceNumber = keyBoardIndices(i);
        break;
    end
end
if deviceNumber == 0
    error('No device by that name was detected');
end

% deviceNumber = 8;

% Setup key press
keyPressNumbers = [KbName('LeftArrow') KbName('RightArrow')];

% Set directories
expDir = pwd; % Set the experimental directory to the current directory 'pwd'
dataDir = 'data'; % Set the path to a directory called 'data'
t.mySeed = sum(100*clock);
rng(t.mySeed); % Make sure to start with a random seed
t.theDate = datestr(now,'yymmdd'); % Collect today's date
t.timeStamp = datestr(now,'HHMM'); % Timestamp

cd(dataDir);
if exist(['vTA_surrSuppressionV2_', p.subject, '.mat'],'file') ~= 0
    load(['vTA_surrSuppressionV2_', p.subject, '.mat']);
    p.runNumber = length(theData)+1;
else
    p.runNumber = 1;
end
cd(expDir);
%% SCREEN PARAMETERS
screens = Screen('Screens'); % look at available screens
p.screenWidthPixels = Screen('Rect', max(screens));
screenWidth = 42; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen, my screen (47.5cm), testing room = 42
viewDistance = 128; % in cm, ideal distance: 1 cm equals 1 visual degree (testing room = 128; my screen = 58)
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % Visual angle of the whole screen
p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
p.grey = 128; white = 255; green = [0 255 0]; blue = [0 0 255]; black = [0 0 0];
%% SOUND SETUP
InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
Fs = 44100;
soundAmp = 1;
reqlatencyclass = 2; % level 2 means take full control over the audio device, even if this causes other sound devices to fail or shutdown
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, Fs, 1); % 1= single-channel
t.cueDur = 0.125; % (s)

cueFreqs = [523.25]; % [high C = precue]
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), t.cueDur, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end
%% GRATING PARAMETERS
p.stimConfigurations = [1 2 3 4 5 6]; 
p.stimConfigurationsNames = {'CollSrrndL' 'CollSrrndR' 'OrthSrrndL' 'OrthSrrndR' 'noSrrndL' 'noSrrndR'};

% contrast parameters
p.fixedStimContrast = 0.2;
p.targContrast = log10(0.6);
p.surroundContrast = 1; 

% size parameters
p.centerSize = round(2 * p.pixPerDeg);
p.surroundSize = round(p.screenWidthPixels(4)); % round(sqrt(p.screenWidthPixels(3)^2 + p.screenWidthPixels(4)^2));
p.gapSize = round(0.02 * p.pixPerDeg);
p.outerFixation = round(0.075 * p.pixPerDeg);
p.ringFixation = p.outerFixation*1.5;
%% TRIAL EVENTS
% create matrix with all unique trial events based on number of repetitions
% 1 = which perception condition - center-surround: center probed or surround probed; center only; surround only; 
% 2 = center contrast
% 3 = surround contrast
% 4 = orientation gratings

[F1] = BalanceFactors(p.numBlocks, 0, p.stimConfigurations);

p.trialEvents = [F1];

p.numTrials = size(p.trialEvents,1);
p.numTrialsPerBlock = p.numTrials/p.numBlocks;
p.numTrialsPerSet = p.numTrialsPerBlock*p.minNumBlocks;
p.numSets = p.numTrials/p.numTrialsPerSet;
p.minNumTrialsPerSC = p.repetitions*2; 

% every trial should be either 45/135; % each target has the same orientation
if mod(p.runNumber,2) ~= 0
    p.targetsOrientation = repmat(45,p.numTrialsPerBlock*p.numBlocks,1);
elseif mod(p.runNumber,2) == 0
    p.targetsOrientation = repmat(135,p.numTrialsPerBlock*p.numBlocks,1); 
end
% surround orientation dependent on collinear/orthogonal condition
p.surroundOrientation = nan(size(p.targetsOrientation));

for nTrial = 1:p.numTrials
   if p.trialEvents(nTrial,1) == 3 || p.trialEvents(nTrial,1) == 4 %  add 90 degrees if orthogonal trial
       p.surroundOrientation(nTrial) = p.targetsOrientation(nTrial) + 90;
   elseif p.trialEvents(nTrial,1) == 1 || p.trialEvents(nTrial,1) == 2 % surround = target if colinear
       p.surroundOrientation(nTrial) = p.targetsOrientation(nTrial);
   end
end

whichOrientation =  [p.targetsOrientation p.surroundOrientation];
p.trialEvents(:, end+1:end+2) = whichOrientation;

% Determine target column
whichTarget = nan(p.numTrials,1);
for nTrial = 1:p.numTrials
    if mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t1 contrast doesn't change (right target)
        whichTarget(nTrial,1) = 2;
    elseif mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t2 contrast doesn't change (left target)
        whichTarget(nTrial,1) = 1;
    end
end
p.trialEvents(:,end+1) = whichTarget;

% Assign cue validity for each trial
p.trialCuesNames = {'Attended' 'Unattended'};

trialCues = zeros(p.numTrials,1);

% assign cues to sets of unique combinations 
for nStimConfig = 1:length(p.stimConfigurations)
   configIndx = find(p.trialEvents(:,1)==p.stimConfigurations(nStimConfig));
   for nAtt = 1:p.numAttTrialsPerComb*p.repetitions
       trialCues(configIndx(nAtt)) = 1;       
   end
   for nUnAtt = 1+(p.numAttTrialsPerComb*p.repetitions):(p.numAttTrialsPerComb*p.repetitions)+((p.minNumBlocks-p.numAttTrialsPerComb)*p.repetitions)
       trialCues(configIndx(nUnAtt)) = 2;
   end
end

attQStructures = 1:9;
unAttQStructures = 10:12;
qStructure = zeros(p.numTrials,1);
p.qStructureNames = {'collCued1' 'orthCued1' 'nsCued1' 'collCued2' 'orthCued2' 'nsCued2' 'collCued3' 'orthCued3' 'nsCued3' 'collUnCued1' 'orthUnCued1' 'nsUnCued1'};
attQCnt = 1;
unAttQCnt = 1;
totAttTrials = p.numAttTrialsPerComb*length(p.stimConfigurations)*p.repetitions;
totUnAttTrials = p.numTrials-totAttTrials;

% assign q structure
for nAtt = 1:totAttTrials
    qStructure(nAtt) = attQStructures(attQCnt);
    if mod(nAtt,2) == 0 && attQCnt < length(attQStructures)
        attQCnt = attQCnt + 1;
    elseif mod(nAtt,2) == 0 && attQCnt >= length(attQStructures)
        attQCnt = 1;
    end
end
for nUnAtt = 1+totAttTrials:totAttTrials+totUnAttTrials
   qStructure(nUnAtt) = unAttQStructures(unAttQCnt);
   if mod(nUnAtt,2) == 0 && unAttQCnt < length(unAttQStructures)
       unAttQCnt = unAttQCnt + 1;
   elseif mod(nUnAtt,2) == 0 && unAttQCnt >= length(unAttQStructures)
       unAttQCnt = 1;
   end
end

p.trialEvents(:,end+1) = qStructure; % store quest structure assignment
p.trialEvents(:,end+1) = trialCues; % store trial cues

% Check trial and cue distribution
trial_cueDistrib = nan(length(p.trialCuesNames), length(p.stimConfigurations)); % [validity x configuration]

for nConfig = 1:length(p.stimConfigurations)
    for nCue = 1:length(p.trialCuesNames)
        if nCue == 1
            trial_cueDistrib(nCue,nConfig) = sum( p.trialEvents(p.trialEvents(:,1)==p.stimConfigurations(nConfig),end-1) == nCue) == p.numAttTrialsPerComb*p.repetitions;
        elseif nCue == 2
            trial_cueDistrib(nCue,nConfig) = sum( p.trialEvents(p.trialEvents(:,1)==p.stimConfigurations(nConfig),end-1) == nCue) == (p.minNumBlocks - p.numAttTrialsPerComb)*p.repetitions;
        end
    end
end

trial_cueDistrib

p.trialEvents % [stimConfiguration, targOrientation, surrOrientation, whichTarget, cueValidity]
p.trialEvents = Shuffle(p.trialEvents,2);
p.trialEvents
p.stimConfigurationsNames

% Define parameters for the stimulus
freq = 2;
p.freq = p.centerSize/p.pixPerDeg * freq;
p.freqSurround = p.surroundSize/p.pixPerDeg * freq;
p.orientation = 0;
numPhases = 2;
p.targPhase = linspace(1,360,numPhases+1);
p.surrPhase = linspace(1,360,numPhases+1);
%% TIMING PARAMETERS
t.startTime = 2; % (s)
t.stimLeadTime = 0.5; % (s)
t.cueTargetSOA = 1; % (s)
t.cueLeadTime = 1; %(s)
t.targetDur = .250; % (s)
t.responseLeadTime = 1; % (s)
t.responseTime = 1; % (s)
t.iti = 0.5; % (s)
% Determine jitter for each trial
t.jit=randsample(0:0.1:1,p.numTrials,true);
t.trialDur = t.jit + (t.cueLeadTime + t.cueTargetSOA + t.targetDur + t.responseLeadTime + t.responseTime + t.iti); % (s)
t.runDur = sum(t.trialDur) + t.startTime*p.numBlocks; % (s)
t.flicker = t.targetDur/2; % (s)
t.numEvents = 7;
t.trialTimes = nan(p.numTrials, t.numEvents); % [startTime cueLeadTime cueTargetSOA targetDur responseLeadTime responseTime iti] [2 1 1 .250 1 1 0.5]
%% CREATE STIMULI
% make mask to create circle for the center grating
[x,y] = meshgrid((-p.centerSize/2):(p.centerSize/2)-1, (-p.centerSize/2):(p.centerSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
centerGaussian = zeros(p.centerSize); centerGaussian(eccen <= (p.centerSize/2)) = 1;
% Gaussian = conv2(Gaussian, fspecial('gaussian', p.pixPerDeg, p.pixPerDeg), 'same');

% Make transparency mask for aplha blending the two images
centerTransparencyMask = zeros(p.centerSize); centerTransparencyMask(eccen <= ((p.centerSize)/2)) = 255;

% make mask to create circle for the surround grating
[x,y] = meshgrid((-p.surroundSize/2):(p.surroundSize/2)-1, (-p.surroundSize/2):(p.surroundSize/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
surroundGaussian = ones(p.surroundSize); % surroundGaussian(eccen <= (p.surroundSize/2)) = 1;

surroundTransparencyMask = zeros(p.surroundSize); surroundTransparencyMask(eccen <= ((p.surroundSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.centerSize-1), 0:(p.centerSize-1));
[Xs,Ys] = meshgrid(0:(p.surroundSize-1), 0:(p.surroundSize-1));

% Make actual gratings
centerGratings = NaN(numPhases, p.centerSize, p.centerSize);
targetGratings = NaN(numPhases,(t.targetDur/t.flicker), p.centerSize, p.centerSize);
surroundGrating = NaN(numPhases, p.surroundSize, p.surroundSize);

for nPhase = 1:numPhases
    center = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.targPhase(nPhase)));
    centerGratings(nPhase,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.surrPhase(nPhase)));
    surroundGrating(nPhase,:,:) = (surround .* surroundGaussian);
end
%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', max(screens), p.grey,[],[],[],[],16);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
load('MyGammaTable.mat');
Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
HideCursor;

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the stimuli
centerX = rect(3)/2; centerY = rect(4)/2;

p.ecc = 3*p.pixPerDeg;

% Coordinates for location on left and right side of fixation
patch =  [centerX-p.ecc centerX+p.ecc]; %[leftStimX rightStimX centerY]
patchDeg = patch/p.pixPerDeg;

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);

t.ifi = Screen('GetFlipInterval',window);
%% Initialize Quest
% Provide our prior knowledge to QuestCreate, and receive the data struct "q".
tGuess = p.targContrast; % intial guess of contrast in log scale
tGuessSd= log10(0.5); % convert to log space!!, Start with a large sd
p.pThreshold = 0.70; % threshold
beta = 2.5; delta = 1/p.numTrials; gamma = 1/2; % gamma 1/amount of answer possibilities

q1 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % collCued1
q2 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % orthCued1
q3 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % nsCued1
q4 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % collCued2
q5 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % orthCued2
q6 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % nsCued2
q7 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % collCued3
q8 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % orthCued3
q9 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % nsCued3
q10 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % collUnCued1
q11= QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % orthUnCued1
q12 = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma); % nsUnCued1

qStructMat = [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12];

q1.normalizePdf = 1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
q2.normalizePdf = 1;
q3.normalizePdf = 1;
q4.normalizePdf = 1;
q5.normalizePdf = 1;
q6.normalizePdf = 1;
q7.normalizePdf = 1;
q8.normalizePdf = 1;
q9.normalizePdf = 1;
q10.normalizePdf = 1;
q11.normalizePdf = 1;
q12.normalizePdf = 1;

%% START THE EXPERIMENT
% Draw some text to the screen first outside of the experimental loop:
% Experiment setup
PsychHID('KbQueueCreate', deviceNumber);
PsychHID('KbQueueStart', deviceNumber);

welcomeText = ['Welcome!' '\n' '\n'...
    '' '\n' '\n'...
    'On every trial two center targets are presented next to fixation, on the left and right.' '\n' '\n' ...
    'Your task is to detect a change in contrast in either one of these center targets. ' '\n' '\n' ...
    'An auditory cue (pre-cue) will indicate when a change will occur.' '\n' '\n' ...
    'On most trials the center targets will be accompanied by a surround stimulus, while on some the surround stimulus is absent.' '\n' '\n' ...
    'Additionally, on some trials there will be no auditory cue.' '\n' '\n' ...
    'When asked to, report which target the change in contrast occured by pressing either the LEFT or RIGHT arrow key.' '\n' '\n' ...
    'Be sure to always maintain steady fixation on the green dot! ' '\n' '\n' '\n' ...
    'Press the DOWN arrow key to continue.' '\n' '\n' ];

DrawFormattedText(window, welcomeText, 'center', 'center', 255);
Screen('Flip', window);
welcomeStart = GetSecs;

% play tones for participant before experiment starts
for nCue=1:size(cueTones,1)
   playSound(pahandle, cueTones(nCue,:)*soundAmp);
   WaitSecs(1);
end

while 1
   [pressed, firstPress] = PsychHID('KbQueueCheck', deviceNumber); % check response
   if pressed == 1
       break;
   end
end

PsychHID('KbQueueStop', deviceNumber);
 
% Preallocate some variables for response data
data.rightwrong = nan(p.numTrials, 1);

data.cThresholdsQ = nan(p.minNumTrialsPerSC,p.numQStructures);  %structure to hold contrast thresholds per trial
qCnt = ones(1,p.numQStructures); % determines which row to access for a structure in cThresholdsQ

% Make surround and center textures before trial loop
surroundStimulus = nan(numPhases,1);
centerStimulus = nan(numPhases,1);
targetStimulus = nan(numPhases,1);

for nPhase = 1:numPhases
    surroundTexture(:,:,1) = squeeze(surroundGrating(nPhase,:,:)) * (p.surroundContrast* p.grey ) + p.grey;
    surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus(nPhase) = Screen('MakeTexture', window, surroundTexture);
    
    centerTexture(:,:,1) = squeeze(centerGratings(nPhase,:,:)) * ( p.fixedStimContrast * p.grey ) + p.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus(nPhase) = Screen('MakeTexture', window, centerTexture);  
end

centerMask = Screen('MakeTexture', window, centerTransparencyMask);

%------------%
% TRIAL LOOP %
%------------%

nSet = 1;
nBlock = 1;

expStart = GetSecs; % baseline experiment start time
welcomeTime = GetSecs - welcomeStart;
t.welcomeTime = welcomeTime;

PsychHID('KbQueueCreate',deviceNumber);
newShift = randsample(1:numPhases,1);

for nTrial = 1:p.numTrials
    currQStruct = p.trialEvents(nTrial,end-1); % which is the currect QUEST structure
    newITI = t.iti + t.jit(nTrial);
    
    % Update contrast threshold with respect to the correct QUEST structure
    p.contrastThreshold = 10^QuestMean(qStructMat(currQStruct));
    data.cThresholdsQ(qCnt(currQStruct),currQStruct) = p.contrastThreshold;
    
    whichTarget = p.trialEvents(nTrial,4);
    
    % cap max contrast
    if p.contrastThreshold > 1
        p.contrastThreshold = 1;
    end    
        p.contrastThreshold
        
    for nPhase = 1:numPhases
        targetTexture(:,:,1) = squeeze(centerGratings(nPhase,:,:)) * ( p.contrastThreshold * p.grey ) + p.grey; 
        targetTexture(:,:,2) = centerTransparencyMask;
        targetStimulus(nPhase) = Screen('MakeTexture', window, targetTexture); 
    end
    
    tic;
    if nTrial == 1 || nTrial == p.numTrialsPerSet*(nSet-1) + 1
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, white, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
        Screen('Flip', window);
        WaitSecs(t.startTime);
    end
    t.trialTimes(nTrial,1) = toc; %duration of start delay
    
    tic;
    nFlicker = 1;
    startCueLeadTime = GetSecs;
    % Turn on stimuli before precue (dur = t.cueLeadTime)
    for nFlips = 1:round(t.cueLeadTime/t.ifi)
        if GetSecs >= (startCueLeadTime + t.cueLeadTime)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
        elseif GetSecs >= startCueLeadTime + t.flicker*nFlicker
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end
        
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, white, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);      
    end
    t.trialTimes(nTrial,2) = toc; % duration of cueLeadTime
    
    % Play pre-cue if valid
    if p.trialEvents(nTrial,end) == 1
        playSound(pahandle, cueTones(1,:)*soundAmp);
        preCueTime = GetSecs - expStart;
    end  
    
    tic;
    nFlicker = 1;
    startCueTargetSOA = GetSecs;
    % Keep stim flickering (dur = t.cueTargetSOA)
    for nFlips = 1:round(t.cueTargetSOA/t.ifi)
        if GetSecs >= (startCueTargetSOA + t.cueTargetSOA)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
        elseif GetSecs >= startCueTargetSOA + t.flicker*nFlicker
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, white, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,3) = toc; %duration of cueTargetSOA
    
    % Set up button press
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush', deviceNumber);
      
    tic;
    nFlicker = 1;
    startTargDur = GetSecs;
    % Flickering Target (dur = t.targetDur)
    for nFlips = 1:round(t.targetDur/t.ifi)
        if GetSecs >= (startTargDur + t.targetDur)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldTargShift = newShift;
            newTargShift = randsample(1:numPhases,1);
            newSurrShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newTargShift == oldTargShift && oldTargShift ~= numPhases
               newTargShift = oldTargShift + 1;
            elseif newTargShift == oldTargShift && oldTargShift == numPhases
                newTargShift = 1;
            end         
        elseif GetSecs >= startTargDur + t.flicker*nFlicker
            oldTargShift = newTargShift;
            newTargShift = randsample(1:numPhases,1);
            oldSurrShift = newSurrShift;
            newSurrShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newTargShift == oldTargShift && oldTargShift ~= numPhases
               newTargShift = oldTargShift + 1;
            elseif newTargShift == oldTargShift && oldTargShift == numPhases
                newTargShift = 1;
            end         
            if newSurrShift == oldSurrShift && oldSurrShift ~= numPhases
               newSurrShift = oldSurrShift + 1;
            elseif newSurrShift == oldSurrShift && oldSurrShift == numPhases
                newSurrShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end
        
        % Draw surroundStimulus if not baseline condition
        if p.trialEvents(nTrial,1) ~= 5 && p.trialEvents(nTrial,1) ~= 6  
            Screen('DrawTexture', window, surroundStimulus(newSurrShift), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], centerX, centerY), p.trialEvents(nTrial,3))
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)')
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)')
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)', p.gapSize, p.gapSize)
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)', p.gapSize, p.gapSize)
        end
        
        % Determine which target to change
        if mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t1 contrast changes (left target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, targetStimulus(newTargShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
            % Draw centerStimulus2
            Screen('DrawTexture', window, centerStimulus(newTargShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2)) 
        elseif mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t2 contrast changes (right target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, centerStimulus(newTargShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
            % Draw centerStimulus2
            Screen('DrawTexture', window, targetStimulus(newTargShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2)) 
        end        
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        % GetClicks;        
        % Stim duration
    end
    t.trialTimes(nTrial,4) = toc; %target duration
    
    tic;
    nFlicker = 1;
    startRespLeadTime = GetSecs;
    % Keep stim flickering (dur = t.responseLeadTime)
    for nFlips = 1:round(t.responseLeadTime/t.ifi)
        if GetSecs >= (startRespLeadTime + t.responseLeadTime)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
        elseif GetSecs >= startRespLeadTime + t.flicker*nFlicker
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,5) = toc; %duration of responseLeadTime
    
    tic;
    nFlicker = 1;
    startRespTime = GetSecs;
    % Keep stim flickering (dur = t.responseTime)
    for nFlips = 1:round(t.responseTime/t.ifi)
        if GetSecs >= (startRespTime + t.responseTime)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
        elseif GetSecs >= startRespTime + t.flicker*nFlicker
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,6) = toc; %duration of responseTime
    
    % Get Behavioral Response
    [pressed, firstPress] = PsychHID('KbQueueCheck', deviceNumber); % check response
    whichPress = find(firstPress);                        
    if any(ismember(whichPress, KbName('ESCAPE')))
        Screen('CloseAll');
        ListenChar(1);
        FlushEvents('keyDown');
        error('User exited program.');                   
    elseif sum(firstPress) == 0 % missed response
        data.rightwrong(nTrial) = 0;
        PsychHID('KbQueueStop', deviceNumber);
    elseif whichPress(1) ~= keyPressNumbers(1) && whichPress(1) ~= keyPressNumbers(2) % irrelevant key press
        data.rightwrong(nTrial) = 0;
        PsychHID('KbQueueStop', deviceNumber);
    elseif (whichPress(1) == keyPressNumbers(1) && whichTarget == 1) || (whichPress(1) == keyPressNumbers(2) && whichTarget == 2) % correct response
        data.rightwrong(nTrial) = 1;
        PsychHID('KbQueueStop', deviceNumber);
    elseif (whichPress(1) == keyPressNumbers(1) && whichTarget == 2) || (whichPress(1) == keyPressNumbers(2) && whichTarget == 1) % incorrect response
        data.rightwrong(nTrial) = 0;
        PsychHID('KbQueueStop', deviceNumber);
    end

    startTrial = GetSecs - expStart; % get the start time of each trial
    
    % Update Quest   
    data.tTestQ(qCnt(currQStruct),currQStruct) = 10^QuestQuantile(qStructMat(currQStruct)); % Recommended by Pelli (1987)
    data.rightwrongQ(qCnt(currQStruct),currQStruct) = data.rightwrong(nTrial);
    qStructMat(currQStruct) = QuestUpdate(qStructMat(currQStruct), log10(data.tTestQ(qCnt(currQStruct),currQStruct)), data.rightwrongQ(qCnt(currQStruct),currQStruct));      
    qCnt(currQStruct) = qCnt(currQStruct) + 1;

    % get ready for next trial   
    tic;
    nFlicker = 1;
    startITI = GetSecs;
    % Keep stim flickering (dur = t.iti)
    for nFlips = 1:round(newITI/t.ifi)
        if GetSecs >= (startITI + newITI)
            break
        end
        % determine new phase if it's time to shift 
        if nFlips == 1
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
        elseif GetSecs >= startITI + t.flicker*nFlicker
            oldShift = newShift;
            newShift = randsample(1:numPhases,1);
            % make sure phase isn't the same as the previous shift
            if newShift == oldShift && oldShift ~= numPhases
               newShift = oldShift + 1;
            elseif newShift == oldShift && oldShift == numPhases
                newShift = 1;
            end
            nFlicker = nFlicker + 1; 
        end    
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(newShift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, black, [centerX-p.ringFixation centerY-p.ringFixation centerX+p.ringFixation centerY+p.ringFixation]);
        Screen('FillOval', window, white, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,7) = toc; % duration of iti
 
    %%% Rest period
    if nTrial == p.numTrialsPerSet*nSet
        rest = GetSecs;       
        restText = ['Set ' num2str(nSet) ' of ' num2str(p.numSets) ' completed! You can take a short break now, ' '' '\n' ...
            'or press the DOWN arrow key to continue' '\n' '\n' ];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);       
        nSet = nSet+1;         
        pmButtonBreak = 0; 
        KbWait; 
        t.restTime(nSet) = (GetSecs-rest)/60;
    end
    
end

t.endTime = GetSecs-expStart; %Get endtime of the experiment in seconds
%Draw some more text to the screen outside of the loop:
Screen(window,'TextSize',30);
byeByeText = 'Great work! You have finished this run.';
DrawFormattedText(window, byeByeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(2);
Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')

% close audio port
PsychPortAudio('Close', pahandle)

% save threshold data
for nQStruct = 1:p.numQStructures
    data.finalThresholdQ(nQStruct) = 10.^QuestMean(qStructMat(nQStruct));
end
%% SAVE OUT THE DATA FILE
cd(dataDir);
theData(p.runNumber).t = t;
theData(p.runNumber).p = p;
theData(p.runNumber).data = data;
% theData(p.runNumber).q = q;
eval(['save vTA_surrSuppression_', p.subject, '.mat theData'])

cd(expDir);