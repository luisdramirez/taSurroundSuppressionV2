
% How does temporal attention (manipulate uncertainty) influence tuned normalization?

%% PREPARE
% Clear workspace, setup trial event parameters (subject name, cue
% validity, repetitions, # of Quest structures), input device, and working
% directories

commandwindow
HideCursor;
echo off
clear all
close all
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 0);

% Subject name
p.subject = 'Pilot_LV';

% Trial Events Parameters
p.cueValidity = 0.75;
[p.numAttTrialsPerComb, p.minNumBlocks] = rat(p.cueValidity);
p.repetitions = 20; %20 reps for at least 40 trials per staircase
p.numBlocks = p.minNumBlocks*p.repetitions; 
p.numQStructures = 12;

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames, ~] = GetKeyboardIndices;
% deviceString = 'Corsair Corsair K95W Gaming Keyboard'; % desk keyboard
deviceString = 'Wired USB Keyboard'; % testing room 207
% deviceString = 'Apple Inc. Apple Keyboard'; % testing room 304
% deviceString = 'Apple Internal Keyboard / Trackpad'; %my laptop
% deviceString = 'CHICONY USB Keyboard'; % yurika's keyboard?

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
% Which screen to display experiment on, pixel size of screen and visual
% angle. Save viewing distance. Save colors. 

screens = Screen('Screens'); % look at available screens
p.screenWidthPixels = Screen('Rect', max(screens));
screenWidth = 34; % my screen =47.5cm; testing room 207 =34; testing room 304 =42 (29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen)
viewDistance = 57; % my screen =58; testing room 207 =57; testing room 304 =128 (in cm, ideal distance: 1 cm equals 1 visual degree) 
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % Visual angle of the whole screen
p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
grey = 128; white = 255; green = [0 255 0]; blue = [0 0 255]; black = [0 0 0]; red = [220 20 60]; 
dimgrey = [105 105 105]; yellow = [255 255 0]; magenta = [255 0 255]; cyan = [0 255 255];

%% GRATING PARAMETERS
% Setup up basic characteristic of gratings: fixed contrast, target
% contrast, surround contrast. Center, surround, gap, and fixation
% size. Stimulus frequency and phase.

% Contrast parameters
p.fixedStimContrast = 0.2; % fixed center contrast 
p.targContrast = log10(0.6); % initial target contrast
p.surroundContrast = 1; % fixed surround contrast

% Size parameters
p.centerSize = round(2 * p.pixPerDeg); % center grating size in pixels
p.surroundSize = round(p.screenWidthPixels(4)); % surround grating size in pixels
p.gapSize = round(0.02 * p.pixPerDeg); % gap between center and surround in pixels
p.fixation = round(0.075 * p.pixPerDeg); % fixation size in pixels
p.fixationRing = p.fixation*1.5; % ring around fixation, dependent on fixation size.

% Frequency and phase parameters
freq = 2;
p.freq = p.centerSize/p.pixPerDeg * freq;
p.freqSurround = p.surroundSize/p.pixPerDeg * freq;
p.orientation = 0;
numPhases = 2;
p.targPhase = linspace(1,360,numPhases+1); % sample random phases for the target gratings
p.surrPhase = linspace(1,360,numPhases+1); % sample random phases for the surround gratings
%% TRIAL EVENTS
% Create matrix with all unique trial events based on number of
% repetitions, saved in p.trialEvents. Stim configurations are the number
% of conditions. This is fed into BalanceFactors. There is only one,
% continuously changing contrast; no need to put it into balance factors.
% Balance factors will output the possible combinations of given
% parameters repeated as many times as wanted. This section then assigns
% target/surround orientation, where the target appears, whether the cue
% will play, and which quest structure a trial will be fed into. Each of
% these is a vector added to p.trialEvents.

%--------------------%
%    Conditions      %
%--------------------%
p.stimConfigurations = [1 2 3 4 5 6]; 
p.stimConfigurationsNames = {'CollSrrndL' 'CollSrrndR' 'OrthSrrndL' 'OrthSrrndR' 'noSrrndL' 'noSrrndR'};
[F1] = BalanceFactors(p.numBlocks, 0, p.stimConfigurations);
p.trialEvents = [F1]; % [stimConfiguration]

p.numTrials = size(p.trialEvents,1);
p.numTrialsPerBlock = p.numTrials/p.numBlocks;
p.numTrialsPerSet = p.numTrialsPerBlock*p.minNumBlocks;
p.numTrialsPerBreak = p.numTrialsPerSet;
p.numSets = p.numTrials/p.numTrialsPerSet;
p.minNumTrialsPerSC = p.repetitions*2; 

%--------------------%
%    Orientations    %
%--------------------%
% Dependent on run number, centers will either by 45/135 deg. 
if mod(p.runNumber,2) ~= 0
    p.targetsOrientation = repmat(45,p.numTrialsPerBlock*p.numBlocks,1);
elseif mod(p.runNumber,2) == 0
    p.targetsOrientation = repmat(135,p.numTrialsPerBlock*p.numBlocks,1); 
end
% Surround orientation is dependent on collinear/orthogonal condition
p.surroundOrientation = nan(size(p.targetsOrientation));
for nTrial = 1:p.numTrials
   if p.trialEvents(nTrial,1) == 3 || p.trialEvents(nTrial,1) == 4 %  add 90 degrees if orthogonal trial
       p.surroundOrientation(nTrial) = p.targetsOrientation(nTrial) + 90;
   elseif p.trialEvents(nTrial,1) == 1 || p.trialEvents(nTrial,1) == 2 % surround = target orientation if colinear
       p.surroundOrientation(nTrial) = p.targetsOrientation(nTrial);
   end
end
whichOrientation =  [p.targetsOrientation p.surroundOrientation]; % Store orientation information
p.trialEvents(:, end+1:end+2) = whichOrientation; % Add orientation info to trialEvents; 
% [stimConfiguration, targOrientation, surrOrientation]

%--------------------%
%      Targets       %
%--------------------%
whichTarget = nan(p.numTrials,1);
for nTrial = 1:p.numTrials
    if mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t1 contrast doesn't change (right target)
        whichTarget(nTrial,1) = 2;
    elseif mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t2 contrast doesn't change (left target)
        whichTarget(nTrial,1) = 1;
    end
end
p.trialEvents(:,end+1) = whichTarget; % Add target info to trialEvents 
% [stimConfiguration, targOrientation, surrOrientation, whichTarget]

%--------------------%
%       Cues         %
%--------------------%
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

%--------------------%
%  Quest Assignment  %
%--------------------%
% Assign Quest structure to each trial
attQStructures = 1:9;
unAttQStructures = 10:12;
qStructure = zeros(p.numTrials,1);
p.qStructureNames = {'collCued1' 'orthCued1' 'nsCued1' 'collCued2' 'orthCued2' 'nsCued2' 'collCued3' 'orthCued3' 'nsCued3' 'collUnCued1' 'orthUnCued1' 'nsUnCued1'};
attQCnt = 1;
unAttQCnt = 1;
totAttTrials = p.numAttTrialsPerComb*length(p.stimConfigurations)*p.repetitions;
totUnAttTrials = p.numTrials-totAttTrials;
% Before shuffling, every 2 trials is 1 condition (1 left and right target trial). These loops assign the
% same structure every 2 trials, depending on whether a trial is an
% attended or unattended trial.
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
p.trialEvents(:,end+1) = qStructure; % add quest structure assignment
p.trialEvents(:,end+1) = trialCues; % add trial cues 
% [stimConfiguration, targOrientation, surrOrientation, whichTarget, qStructures, cueValidity]

%---------------------%
% Check Distributions %
%---------------------%
% trialCueDistrib should be all 1s (num of attended and unattended should be equal to the cue validity ratio*numRepetions). 
trial_cueDistrib = nan(length(p.trialCuesNames), length(p.stimConfigurations)); % [validity x configuration]
for nConfig = 1:length(p.stimConfigurations)
    for nCue = 1:length(p.trialCuesNames)
        if nCue == 1
            trial_cueDistrib(nCue,nConfig) = sum( p.trialEvents(p.trialEvents(:,1)==p.stimConfigurations(nConfig),end) == nCue) == p.numAttTrialsPerComb*p.repetitions;
        elseif nCue == 2
            trial_cueDistrib(nCue,nConfig) = sum( p.trialEvents(p.trialEvents(:,1)==p.stimConfigurations(nConfig),end) == nCue) == (p.minNumBlocks - p.numAttTrialsPerComb)*p.repetitions;
        end
    end
end
trial_cueDistrib; 
p.trialEvents = Shuffle(p.trialEvents,2);
p.trialEvents; % [stimConfiguration, targOrientation, surrOrientation, whichTarget, qStructures, cueValidity]
p.stimConfigurationsNames;
%% TIMING PARAMETERS

% Setup basic timing of events.
t.startTime = 2; % (s) delay before stimulus comes up when a new set starts (after a break)
t.cueLeadTime = .250; %(s) time between onset of stimulus and cue
t.cueTargetSOA = .500; % (s) time between cue and target
t.targetDur = .250; % (s) target duration
t.responseTime = 1; % (s) response duration after target appears
t.feedback = 0.250; % (s) feedback duration
t.iti = 0.0; % (s) dwell period between trials

% Determine jitter for each trial
t.jit=randsample([0:0.25:2],p.numTrials,true); % jitter to add to iti; ranges from nothing to 1s
t.trialDur = t.jit + (t.cueLeadTime + t.cueTargetSOA + t.targetDur + t.responseTime + t.iti); % (s) each trial duration
t.runDur = sum(t.trialDur) + t.startTime*p.numBlocks; % (s) total run duration (does not include rest time)
t.flicker = t.targetDur/2; % (s) how many times to shift phase.

% Allocate space to record timing of events.
t.numEvents = 7;
t.trialTimes = nan(p.numTrials, t.numEvents); % [startTime cueLeadTime cueTargetSOA targetDur responseLeadTime responseTime iti] [2 1 1 .250 1 1 0.5]
%% CREATE STIMULI
% Generate center and surround grating, in addition to their masks. 

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
% Open window and where to display center targets.

[window,rect] = Screen('OpenWindow', max(screens), grey,[],[],[],[],16);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
load('linearizedCLUT.mat');
Screen('LoadNormalizedGammaTable', window, linearizedCLUT);
HideCursor;

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the stimuli
centerX = rect(3)/2; centerY = rect(4)/2; % center coordinates
p.ecc = 3*p.pixPerDeg; % eccentricity of targets from center 

% Coordinates for location on left and right side of fixation
patch =  [centerX-p.ecc centerX+p.ecc]; %[leftStimX rightStimX] % X coordinates of targets
patchDeg = patch/p.pixPerDeg; % coordinates of targets in degrees

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
t.ifi = Screen('GetFlipInterval',window); % grab screen refresh rate
%% SOUND SETUP
% Initialize audio channel to play cues. Sets up cues by number of elements
% in cueFreqs, which determines which frequency to play.

InitializePsychSound(1); % 1 for precise timing

% Open audio device for low-latency output
Fs = 44100;
soundAmp = 1;
reqlatencyclass = 2; % level 2 means take full control over the audio device, even if this causes other sound devices to fail or shutdown
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, Fs, 1); % 1= single-channel

cueFreqs = 10.^linspace(log10(261.63),log10(880),round(t.cueTargetSOA/t.ifi)); % 
%Create tones
for iTone = 1:numel(cueFreqs)
    tone = MakeBeep(cueFreqs(iTone), t.ifi, Fs);
    cueTones(iTone,:) = applyEnvelope(tone, Fs);
end
%% INITIALIZE QUEST
% Give intial contrast, threshold to begin Quest
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

qStructVect = [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12];

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

data.qMissedTrials = zeros(1,p.numQStructures); % vector to keep track of missed trials in a quest structure
data.qBadPresses = zeros(1,p.numQStructures); % vector to keep track of bad button presses in a quest structure
%% EXPERIMENT 
% Draw some text to the screen first outside of the experimental loop:
% Experiment setup
PsychHID('KbQueueCreate', deviceNumber);
PsychHID('KbQueueStart', deviceNumber);

welcomeText = ['Welcome!' '\n' '\n'...
    '' '\n' '\n'...
    'On every trial two center targets are presented next to fixation, on the left and right.' '\n' '\n' ...
    'Your task is to detect a change in contrast in either one of these center targets. ' '\n' '\n' ...
    'An auditory cue (pre-cue) will indicate when an increment in contrast will occur.' '\n' '\n' ...
    'On most trials the center targets will be accompanied by a surround stimulus, while on some the surround stimulus is absent.' '\n' '\n' ...
    'Additionally, on some trials there will be no auditory cue.' '\n' '\n' ...
    'When the fixation becomes green, report which target the contrast increment occured by pressing either the LEFT or RIGHT arrow key.' '\n' '\n' ...
    'Be sure to always maintain steady fixation on dot! ' '\n' '\n' '\n' ...
    'Press the DOWN arrow key to continue.' '\n' '\n' ];

DrawFormattedText(window, welcomeText, 'center', 'center', 255);
Screen('Flip', window);
welcomeStart = GetSecs;

% Play tones for participant before experiment starts
for nCue=1:size(cueTones,1)
   playSound(pahandle, cueTones(nCue,:)*soundAmp);
end

while 1
   [pressed, firstPress] = PsychHID('KbQueueCheck', deviceNumber); % check response
   if pressed == 1
       break;
   end
end

PsychHID('KbQueueStop', deviceNumber);
 
% Preallocate some variables to store response data
data.rightwrong = nan(p.numTrials, 1); % collect behavioral response
data.cThresholdsQ = nan(p.minNumTrialsPerSC,p.numQStructures);  %structure to hold contrast thresholds per trial
qCnt = ones(1,p.numQStructures); % determines which row to access for a structure in cThresholdsQ
missed = zeros(p.numTrials,1); badPress = zeros(p.numTrials,1);
 
% Make surround and center textures before trial loop
surroundStimulus = nan(numPhases,1);
centerStimulus = nan(numPhases,1);
targetStimulus = nan(numPhases,1);

for nPhase = 1:numPhases
    surroundTexture(:,:,1) = squeeze(surroundGrating(nPhase,:,:)) * (p.surroundContrast* grey ) + grey;
    surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus(nPhase) = Screen('MakeTexture', window, surroundTexture);
    
    centerTexture(:,:,1) = squeeze(centerGratings(nPhase,:,:)) * ( p.fixedStimContrast * grey ) + grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus(nPhase) = Screen('MakeTexture', window, centerTexture);  
end

centerMask = Screen('MakeTexture', window, centerTransparencyMask);

%------------%
% TRIAL LOOP %
%------------%
% Each loop will break once the desired time for a section of an experiment
% is met. A flip is dependent on the refresh rate in a given window of
% time. A flicker (phase shift update) occurs when the duration of a
% flicker passes. Each loop ensures each phase is different from the
% previous. 

% Initialize timing, counters, some vars
nSet = 0; % set counter
nBreak = 1; %rest counter
expStart = GetSecs; % baseline experiment start time
welcomeTime = GetSecs - welcomeStart;
t.welcomeTime = welcomeTime;
nTrial = 1; % completed trials, only updated if previous trial is successful
currTrial = 0; % attempted trials, updates every loop
newShift = randsample(1:numPhases,1); % prelim phase shift

PsychHID('KbQueueCreate',deviceNumber);

while nTrial <= p.numTrials
    
    % Check if retry is necessary
    if missed(nTrial) == 1 || badPress(nTrial) == 1
        missed(nTrial) = 0; badPress(nTrial) = 0;
    end
    
    %--------------------%
    %    Index Updates   %
    %--------------------%    
    whichTarget = p.trialEvents(nTrial,4);
    whichQStruct = p.trialEvents(nTrial,end-1); % which is the currect QUEST structure
    newITI = t.iti + t.jit(nTrial); % what is this trial's new iti

    %--------------------%
    %  Threshold Update  %
    %--------------------%
    % Update contrast threshold with respect to the correct QUEST structure
    p.contrastThreshold = 10^QuestMean(qStructVect(whichQStruct));
    data.cThresholdsQ(qCnt(whichQStruct),whichQStruct) = p.contrastThreshold;
    % Cap max contrast
    if p.contrastThreshold > 1
        p.contrastThreshold = 1;
    end    
    p.contrastThreshold;
    
    %--------------------%
    %    Target Update   %
    %--------------------%
    for nPhase = 1:numPhases
        targetTexture(:,:,1) = squeeze(centerGratings(nPhase,:,:)) * ( p.contrastThreshold * grey ) + grey; 
        targetTexture(:,:,2) = centerTransparencyMask;
        targetStimulus(nPhase) = Screen('MakeTexture', window, targetTexture); 
    end
    
    %--------------------%
    %      Begin Exp     %
    %--------------------%
    % Draw just the fixation at the beginning, or after a break
    tic;
    if currTrial == 0 || (currTrial > 1 && currTrial == p.numTrialsPerBreak*(nBreak-1))
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, white, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation]);
        Screen('Flip', window);
        WaitSecs(t.startTime);
    end
    t.trialTimes(nTrial,1) = toc; % duration of start delay
    
    %--------------------%
    %      Before Cue    %
    %--------------------%
    % Turn on stimuli before precue (dur = t.cueLeadTime)
    tic;
    nFlicker = 1;
    startCueLeadTime = GetSecs;
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, white, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);      
    end
    t.trialTimes(nTrial,2) = toc; % duration of cueLeadTime
      
    %--------------------%
    %         Cue        % 
    %--------------------%
    % Keep stim flickering after cue plays (dur = t.cueTargetSOA)   
    tic;
    nFlicker = 1;
    startCueTargetSOA = GetSecs;
    for nFlips = 1:round(t.cueTargetSOA/t.ifi)
        if GetSecs >= (startCueTargetSOA + t.cueTargetSOA)
            break
        end
        % Play pre-cue if valid!
        if p.trialEvents(nTrial,end) == 1
            playSound(pahandle, cueTones(nFlips,:)*soundAmp);          
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, white, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,3) = toc; % duration of cueTargetSOA
    
    % Prepare button press for behavioral response
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush', deviceNumber);
        
    %--------------------%
    %       Target       %
    %--------------------%
    % Display target, keeping stimuli flickering (dur = t.targetDur)  
    tic;
    nFlicker = 1;
    startTargDur = GetSecs;
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
            Screen('FillOval', window, grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)')
            Screen('FillOval', window, grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)')
            Screen('FrameOval', window, grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)', p.gapSize, p.gapSize)
            Screen('FrameOval', window, grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)', p.gapSize, p.gapSize)
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, black, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);
        % GetClicks;        
        % Stim duration
    end
    t.trialTimes(nTrial,4) = toc; % duration of target
    
    %--------------------%
    %      Response      %
    %--------------------%
    % Keep stim flickering for response time, fixation becomes green (dur = t.responseTime)    
    tic;
    nFlicker = 1;
    startRespTime = GetSecs;
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, black, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,5) = toc; % duration of responseTime
    
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
        missed(nTrial) = 1;
        data.qMissedTrials(whichQStruct) = data.qMissedTrials(whichQStruct) + 1;
        PsychHID('KbQueueStop', deviceNumber);
    elseif whichPress(1) ~= keyPressNumbers(1) && whichPress(1) ~= keyPressNumbers(2) % irrelevant key press
        data.rightwrong(nTrial) = 0;
        badPress(nTrial) = 1;
        data.qBadPresses(whichQStruct) = data.qBadPresses(whichQStruct) + 1;
        PsychHID('KbQueueStop', deviceNumber);
    elseif (whichPress(1) == keyPressNumbers(1) && whichTarget == 1) || (whichPress(1) == keyPressNumbers(2) && whichTarget == 2) % correct response
        data.rightwrong(nTrial) = 1;
        PsychHID('KbQueueStop', deviceNumber);
    elseif (whichPress(1) == keyPressNumbers(1) && whichTarget == 2) || (whichPress(1) == keyPressNumbers(2) && whichTarget == 1) % incorrect response
        data.rightwrong(nTrial) = 0;
        PsychHID('KbQueueStop', deviceNumber);
    end
    
    %--------------------%
    %    Update Quest    %
    %--------------------%
    if missed(nTrial) == 0 && badPress(nTrial) == 0 % if the trial was successful, update quest
        data.tTestQ(qCnt(whichQStruct),whichQStruct) = 10^QuestQuantile(qStructVect(whichQStruct)); % Recommended by Pelli (1987)
        % insert repsonse to current quest structure
        data.rightwrongQ(qCnt(whichQStruct),whichQStruct) = data.rightwrong(nTrial); 
        % update current quest structure 
        qStructVect(whichQStruct) = QuestUpdate(qStructVect(whichQStruct), log10(data.tTestQ(qCnt(whichQStruct),whichQStruct)), data.rightwrongQ(qCnt(whichQStruct),whichQStruct));       
        % update the number of times the current quest structure has been updated (keeps track of how many trials are in the structure)
        qCnt(whichQStruct) = qCnt(whichQStruct) + 1; 
    end
    
    %--------------------%
    %      Feedback      %
    %--------------------%
    % Determine color of fixation feedback
    if data.rightwrong(nTrial) == 1 % correct = green
        feedbackColor = green;
    elseif data.rightwrong(nTrial) == 0 && (missed(nTrial) == 0 && badPress(nTrial) == 0) == 1 % incorrect = red
        feedbackColor = red;
    elseif missed(nTrial) == 1 % missed = blue
        feedbackColor = blue;
    elseif badPress(nTrial) == 1 % bad press = black
        feedbackColor = grey;
    end
    % Keep stim flickering for feedback time, fixation becomes green (dur = t.feedback)
    tic;
    nFlicker = 1;
    startFeedback = GetSecs;
    for nFlips = 1:round(t.feedback/t.ifi)
        if GetSecs >= (startFeedback + t.feedback)
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
        elseif GetSecs >= startFeedback + t.flicker*nFlicker
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, feedbackColor, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,6) = toc; % duration of feedback
    
    %--------------------%
    %        ITI         %
    %--------------------%   
    % Get ready for next trial, keep stim flickering (dur = newITI)
    tic;
    nFlicker = 1;
    startITI = GetSecs;
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
        Screen('FillOval', window, black, [centerX-p.fixationRing centerY-p.fixationRing centerX+p.fixationRing centerY+p.fixationRing]);
        Screen('FillOval', window, white, [centerX-p.fixation centerY-p.fixation centerX+p.fixation centerY+p.fixation])
        Screen('Flip', window);
    end
    t.trialTimes(nTrial,7) = toc; % duration of iti
    
    %--------------------%
    %       Updates      %   
    %--------------------%
    % if no miss or bad press, update trial index || reshuffle trial into trial order 
    if missed(nTrial) == 0 && badPress(nTrial) == 0
        nTrial = nTrial+1;
    elseif missed(nTrial) == 1 || badPress(nTrial) == 1
        currTrialInfo = p.trialEvents(nTrial,:); % store current trial
        p.trialEvents(nTrial,:) = []; % delete current trial from trial events
        p.trialEvents(end+1,:) = currTrialInfo; % append current trial to end of trial events
        p.trialEvents = [p.trialEvents(1:nTrial-1,:); Shuffle(p.trialEvents(nTrial:end,:),2)]; % reshuffle trial events
        p.trialEvents; % verify shuffle
    end    
    % Update number of attempted trials
    currTrial = currTrial+1;   
    % Update sets completed
    if mod(nTrial,p.numTrialsPerSet) == 0
        nSet = nSet+1;
    end
    
    %--------------------%
    %     Rest Period    %
    %--------------------%
    if mod(currTrial,p.numTrialsPerBreak) == 0       
        rest = GetSecs;       
        restText = ['You can take a short break now, ' '' '\n' ...
            'or press the DOWN arrow key to continue' '\n' '\n'...
            'Progress: set ' num2str(nSet) ' of ' num2str(p.numSets) ' complete.'];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);                
        pmButtonBreak = 0; 
        KbWait; 
        t.restTime(nBreak) = (GetSecs-rest)/60;
        nBreak = nBreak+1;
    end
    
end

%--------------------%
%   End Experiment   %
%--------------------%
t.endTime = GetSecs-expStart; % Get endtime of the experiment in seconds

% Draw some more text to the screen outside of the loop, close experiment
Screen(window,'TextSize',30);
byeByeText = 'Great work! You have finished this run.';
DrawFormattedText(window, byeByeText, 'center', 'center', [255 255 255]);
Screen('Flip', window);
WaitSecs(2);
Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
Screen('CloseAll')

% Close audio port
PsychPortAudio('Close', pahandle)

% Save threshold data
for nQStruct = 1:p.numQStructures
    data.finalThresholdQ(nQStruct) = 10.^QuestMean(qStructVect(nQStruct));
end
%% SAVE OUT THE DATA FILE
cd(dataDir);
theData(p.runNumber).t = t;
theData(p.runNumber).p = p;
theData(p.runNumber).data = data;
eval(['save vTA_surrSuppressionV2_', p.subject, '.mat theData'])
cd(expDir);