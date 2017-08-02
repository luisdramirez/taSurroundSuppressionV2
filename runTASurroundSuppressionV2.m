
% How does alertness/arousal/temporal attention influence tuned
% normalization cruves?

% If attention enhances orientation-tuned suppression, measures of tuned
% normalization bandwidth may become narrower when attention is directed
% towards a stimulus.

% 

%% PREPARE AND COLLECT INFO

commandwindow
HideCursor;
echo off
clear all
close all
KbName('UnifyKeyNames');
Screen('Preference', 'SkipSyncTests', 0);

% Subject name and run number

p.subject = 'Pre-Pilot_LR';
p.cueValidity = 0.75;
[p.numAttTrialsPerComb, p.minNumBlocks] = rat(p.cueValidity);
p.repetitions = 8;
p.numBlocks = p.minNumBlocks*p.repetitions; 
p.numBreaks = p.numBlocks*2;

useEyeTracker = 'No';

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
deviceNumber = 8;
% Setup key presse
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
screenWidth = 42; % 29 cm macbook air, 40 cm trinitron crt, 60 cm Qnix screen
viewDistance = 128; % in cm, ideal distance: 1 cm equals 1 visual degree
visAngle = (2*atan2(screenWidth/2, viewDistance))*(180/pi); % Visual angle of the whole screen
p.pixPerDeg = round(p.screenWidthPixels(3)/visAngle); % pixels per degree visual angle
p.grey = 128;
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

%playSound(pahandle, cueTones(1,:)*soundAmp);
%playSound(pahandle, cueTones(2,:)*soundAmp);


%% GRATING PARAMETERS
p.stimConfigurations = [1 2 3 4 5 6]; 
p.stimConfigurationsNames = {'CollSrrndL' 'CollSrrndR' 'OrthSrrndL' 'OrthSrrndR' 'noSrrndL' 'noSrrndR'};

% contrast parameters
p.fixedStimContrast = 0.2;
p.targContrast = log10(0.75);
p.surroundContrast = 1; 

% size parameters
p.centerSize = round(1 * p.pixPerDeg);
p.surroundSize = round(sqrt(p.screenWidthPixels(3)^2 + p.screenWidthPixels(4)^2));
p.gapSize = round(0.02 * p.pixPerDeg);
p.outerFixation = round(0.05 * p.pixPerDeg);
p.innerFixation = p.outerFixation/1.5;


%% TRIAL EVENTS
% create matrix with all unique trial events based on number of repetitions
% 1 = which perception condition - center-surround: center probed or surround probed; center only; surround only; 
% 2 = center contrast
% 3 = surround contrast
% 4 = orientation gratings

[F1] = BalanceFactors(p.numBlocks, 0, p.stimConfigurations);

p.trialEvents = [F1];

p.numTrialsPerConfig = sum(p.trialEvents(:,1) == 1);

p.numTrials = size(p.trialEvents,1);
p.numTrialsPerBreak = p.numTrials/p.numBreaks;
p.numTrialsPerBlock = p.numTrials/p.numBlocks;

% every trial should be a either 45/135; 
p.targetsOrientation = [repmat(45,p.numTrialsPerBlock*(p.numBlocks/2),1);repmat(135,p.numTrialsPerBlock*(p.numBlocks/2),1)]; % each target has the same orientation

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

p.trialEvents(:,end+1) = trialCues; % store trial cues at the end trialEvents

% Check trial and cue distribution
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
numPhases = 20;
p.phase = randsample(1:180,numPhases*2, true); % [centerPhase surrPhase]
p.phase = reshape(p.phase, [numPhases 2]);
%% TIMING PARAMETERS
t.targetDur = .250; % (s)
t.iti = 1; % (s)
t.startTime = 2; % (s)
t.stimLeadTime = 1; % (s)
t.responseTime = []; % (s)
t.cueTargetSOA = 1; % (s)
t.cueLeadTime = 1; %(s)
t.responseLeadTime = 1; % (s)
t.trialDur = t.cueLeadTime + t.cueTargetSOA*2 + t.targetDur*2 + t.responseLeadTime; % duration of the longest trial (sS)
t.trialDurLongest = t.trialDur + t.startTime; % (2)
t.runDur = t.trialDur*p.numTrials; % (s)
t.flicker = 0.025; % (s)
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

% surroundTransparencyMask = zeros(p.surroundSize); surroundTransparencyMask(eccen >= ((p.centerSize)/2)) = 255;

% make unique grating for every trial
[Xc,Yc] = meshgrid(0:(p.centerSize-1), 0:(p.centerSize-1));
[Xs,Ys] = meshgrid(0:(p.surroundSize-1), 0:(p.surroundSize-1));

% Make actual gratings
centerGratings = NaN(p.numTrials*2, p.centerSize, p.centerSize);
targetGratings = NaN(p.numTrials*2,(t.targetDur/t.flicker), p.centerSize, p.centerSize);
surroundGrating = NaN(p.numTrials, p.surroundSize, p.surroundSize);

for nPhase = 1:numPhases
    center = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(nPhase,1)));
    centerGratings(nPhase,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(nPhase,2)));
    surroundGrating(nPhase,:,:) = (surround .* surroundGaussian);
end

%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', max(screens), p.grey,[],[],[],[],16);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
load('MyGammaTable.mat');
Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
HideCursor;
white = 255; green = [0 255 0]; blue = [0 0 255];

% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the stimuli
centerX = rect(3)/2; centerY = rect(4)/2;
% Coordinates for location on left and right side of fixation
patch =  [centerX*(4/5) centerX*(6/5)]; %[leftStimX rightStimX centerY]

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);

%% EYE TRACKING
if strcmp(useEyeTracker, 'Yes')
    p.observer = [p.subject num2str(p.runNumber)];
    [el edf_filename] = eyeTrackingOn(window, p.observer, rect, p.pixPerDeg);
end

% [status] = Eyelink('Message', string, trigger) % prints out message
% [status] = Eyelink('CheckRecording') % 0 if recording, 1 if not recording 
% [time] = Eyelink('TrackerTime') % time since the tracker application started
%% START THE EXPERIMENT

% Initialize Quest

% Provide our prior knowledge to QuestCreate, and receive the data struct "q".
tGuess = p.targContrast; % intial guess of contrast in log scale
tGuessSd= log10(0.5); % convert to log space!!, Start with a large sd
p.pThreshold = 0.80; % threshold
beta = 2.5; delta = 1/p.numTrials; gamma = 1/2; % gamma 1/amount of answer possibilities
qSurr = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma);
qSurr.normalizePdf = 1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
qNSurr = QuestCreate(tGuess, tGuessSd, p.pThreshold, beta, delta, gamma);
qNSurr.normalizePdf = 1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

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
   if pressed == 1;
       break;
   end
end

PsychHID('KbQueueStop', deviceNumber);
 
% Preallocate some variables for response data
data.rightwrong = zeros(p.numTrials,1);
data.cThresholdsSurr = nan(p.numTrials,1);
data.cThresholdsNSurr = nan(p.numTrials,1);

% Make surround and center textures before trial loop
surroundStimulus = nan(numPhases,1);
centerStimulus = nan(numPhases,1);
targetStimulus = nan(numPhases,1);

for nPhase = 1:numPhases
    surroundTexture(:,:,1) = squeeze(surroundGrating(nPhase,:,:)) * (p.surroundContrast* p.grey ) + p.grey;
%     surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus(nPhase) = Screen('MakeTexture', window, surroundTexture);
    
    centerTexture(:,:,1) = squeeze(centerGratings(nPhase,:,:)) * ( p.fixedStimContrast * p.grey ) + p.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus(nPhase) = Screen('MakeTexture', window, centerTexture);  
end

centerMask = Screen('MakeTexture', window, centerTransparencyMask);

%------------%
% TRIAL LOOP %
%------------%

% Eye Tracker Recording ON
if strcmp(useEyeTracker, 'Yes')
    eyeTrackingRecord(el, rect, p.pixPerDeg);
end

nBlock = 1;
nBreak = 1;

expStart = GetSecs; % baseline experiment start time
welcomeTime = GetSecs - welcomeStart;
t.welcomeTime = welcomeTime;

PsychHID('KbQueueCreate',deviceNumber);

for nTrial = 1:p.numTrials
    
    if p.trialEvents(nTrial,1) ~= 5 && p.trialEvents(nTrial,1) ~= 6 %use surround quest structure
        p.contrastThreshold = 10^QuestMean(qSurr);
        data.cThresholdsSurr(nTrial) = p.contrastThreshold;
    elseif p.trialEvents(nTrial,1) == 5 || p.trialEvents(nTrial,1) == 6
        p.contrastThreshold = 10^QuestMean(qNSurr); % use no surround ques structure
        data.cThresholdsNSurr(nTrial) = p.contrastThreshold;
    end
    
    whichTarget = p.trialEvents(nTrial, end-1);
    
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
    
    if nTrial == 1 || nTrial == p.numTrialsPerBlock*(nBlock-1) + 1
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation]);
        Screen('Flip', window);
        WaitSecs(t.startTime);
    end
    
    % Draw Fixation
    Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    WaitSecs(t.stimLeadTime);
    
    % Turn on stimuli before precue (dur = t.cueLeadTime)
    for nShift = 1:(t.cueLeadTime/t.flicker)       
        shift = randsample(1:numPhases,1);
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);       
        WaitSecs(t.flicker);
    end
    
    % Play pre-cue if valid
    if p.trialEvents(nTrial,end) == 1
        playSound(pahandle, cueTones(1,:)*soundAmp);
        preCueTime = GetSecs - expStart;
        trialTimes(nTrial,2) = preCueTime;
    end  
    
    % Keep stim flickering (dur = t.cueTargetSOA)
    for nShift = 1:(t.cueTargetSOA/t.flicker)  
        shift = randsample(1:numPhases,1);      
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end
    
    % Set up button press
    PsychHID('KbQueueStart', deviceNumber);
    PsychHID('KbQueueFlush', deviceNumber);
    
    % Flickering Target (dur = t.targetDur)
    for nShift = 1:(t.targetDur/t.flicker)
        shift = randsample(1:numPhases,1);
        % Draw surroundStimulus if not baseline condition
        if p.trialEvents(nTrial,1) ~= 5 && p.trialEvents(nTrial,1) ~= 6  
            Screen('DrawTexture', window, surroundStimulus(shift), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], centerX, centerY), p.trialEvents(nTrial,3))
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)')
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)')
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)', p.gapSize, p.gapSize)
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)', p.gapSize, p.gapSize)
        end
        
        % Determine which target to change
        if mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t1 contrast changes (left target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, targetStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
            % Draw centerStimulus2
            Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2)) 
        elseif mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t2 contrast changes (right target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
            % Draw centerStimulus2
            Screen('DrawTexture', window, targetStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2)) 
        end
        
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        % GetClicks;        
        % Stim duration
        WaitSecs(t.flicker);
    end

    % Keep stim flickering (dur = t.responseLeadTime)
    for nShift = 1:(t.responseLeadTime/t.flicker)
        shift = randsample(1:numPhases,1);
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,2))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,2))
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end

    Screen('FillOval', window, blue, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    
    % Get Behavioral Response
    check = 0;
    while 1
        [pressed, firstPress] = PsychHID('KbQueueCheck', deviceNumber); % check response
        if pressed
            whichPress = find(firstPress);
            
            if check == 0
                
                if any(ismember(whichPress, KbName('ESCAPE')))
                    Screen('CloseAll');
                    ListenChar(1);
                    FlushEvents('keyDown');
                    error('User exited program.');                   
                elseif (whichPress(1) == keyPressNumbers(1) && whichTarget == 1) || (whichPress(1) == keyPressNumbers(2) && whichTarget == 2) % correct response
                    data.rightwrong(nTrial) = 1;
                    check = 1;
                    PsychHID('KbQueueStop', deviceNumber);
                else
                    check = 1;
                    PsychHID('KbQueueStop', deviceNumber);
                end
            end
            %exit loop if keypress is detected
            if check == 1
                break;
            end
        end
        

    end
    startTrial = GetSecs - expStart; % get the start time of each trial
    
    % Update Quest
    if p.trialEvents(nTrial,1) ~= 5 && p.trialEvents(nTrial,1) ~= 6 % use surr quest structure
        data.tTest(nTrial) = 10^QuestQuantile(qSurr); % Recommended by Pelli (1987)
        qSurr = QuestUpdate(qSurr, log10(data.tTest(nTrial)), data.rightwrong(nTrial));
    elseif p.trialEvents(nTrial,1) == 5 || p.trialEvents(nTrial,1) == 6 % use no surr quest structure
        data.tTest(nTrial) = 10^QuestQuantile(qNSurr); % Recommended by Pelli (1987)
        qNSurr = QuestUpdate(qNSurr, log10(data.tTest(nTrial)), data.rightwrong(nTrial));        
    end
    % get ready for next trial
    Screen('Flip', window);
    WaitSecs(t.iti);
 
    %%% Rest period
    if nTrial == p.numTrialsPerBreak*nBreak && nTrial ~= p.numTrialsPerBlock*nBlock
        rest = GetSecs;       
        restText = 'You can take a short break now, or press the DOWN arrow key to continue.';
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);       
        nBreak = nBreak+1;        
        pmButtonBreak = 0;         
        KbWait;        
        t.restTime = (GetSecs-rest)/60;
    elseif nTrial == p.numTrialsPerBlock*nBlock
        rest = GetSecs;       
        restText = ['Block ' num2str(nBlock) ' of ' num2str(p.numBlocks) ' completed! You can take a short break now, ' '' '\n' ...
            'or press the DOWN arrow key to continue' '\n' '\n' ];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);       
        nBlock = nBlock+1;
        nBreak = nBreak+1;         
        pmButtonBreak = 0; 
        KbWait; 
        t.restTime = (GetSecs-rest)/60;
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

% Eye Tracker Recording OFF
if strcmp(useEyeTracker, 'Yes')
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    Eyelink('ReceiveFile',edf_filename);
end

% close audio port
PsychPortAudio('Close', pahandle)
data.finalThresholdSurr = 10.^QuestMean(qSurr);
data.finalThresholdNSurr = 10.^QuestMean(qNSurr);

data.cThresholdsSurr = data.cThresholdsSurr(~isnan(data.cThresholdsSurr));
data.cThresholdsNSurr = data.cThresholdsNSurr(~isnan(data.cThresholdsNSurr));

%% SAVE OUT THE DATA FILE
cd(dataDir);
theData(p.runNumber).t = t;
theData(p.runNumber).p = p;
theData(p.runNumber).data = data;
% theData(p.runNumber).q = q;
eval(['save vTA_surrSuppression_', p.subject, '.mat theData'])

cd(expDir);
