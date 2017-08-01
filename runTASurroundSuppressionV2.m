
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
p.repetitions = 1;
p.numBlocks = p.minNumBlocks*p.repetitions; 
p.numBreaks = p.numBlocks*2;

useEyeTracker = 'No';

% Check which devicenumber the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames] = GetKeyboardIndices;
deviceString = 'Corsair Corsair K95W Gaming Keyboard';
% deviceString = 'Apple Inc. Apple Keyboard';
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
p.numContrasts = 1;
p.minContrast = 0.1;
p.maxContrast = 0.75;
p.fixedStimContrast = 0.2; 
p.targContrast = 0.75; % the increment contrast
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

[F1, F2, F3] = BalanceFactors(p.numBlocks, 0, p.stimConfigurations, p.targContrast, p.targContrast);

p.trialEvents = [F1, F2, F3];

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

% Determine target contrast
for nTrial = 1:p.numTrials
    if mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t1 contrast doesn't change (right target)
        p.trialEvents(nTrial,2) = p.fixedStimContrast;
    elseif mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t2 contrast doesn't change (left target)
        p.trialEvents(nTrial,3) = p.fixedStimContrast;
    end
end

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

p.trialEvents % [stimConfiguration, t1Contrast, t2Contrast, targOrientation, surrOrientation, cueValidity]
% p.trialEvents = Shuffle(p.trialEvents,2);
p.trialEvents
p.stimConfigurationsNames 

% Define parameters for the stimulus
freq = 2;
p.freq = p.centerSize/p.pixPerDeg * freq;
p.freqSurround = p.surroundSize/p.pixPerDeg * freq;
p.orientation = 0;
p.phase = randsample(1:180,p.numTrials*2, true); % [centerPhase surrPhase]
p.phase = reshape(p.phase, [p.numTrials 2]);
%% TIMING PARAMETERS
t.targetDur = .5; % (s)
t.iti = 1; % (s)
t.startTime = 2; % (s)
t.stimLeadTime = 1; % (s)
t.responseTime = []; % (s)
t.cueTargetSOA = 1; % (s)
t.cueLeadTime = 1; %(s)
t.responseLeadTime = 1; % (s)
t.trialDur = t.cueLeadTime + t.cueTargetSOA*2 + t.targetDur*2 + t.responseLeadTime; % duration of the longest trial
t.trialDurLongest = t.trialDur + t.startTime;

jit = 0:0.2:1;
trialJit = Shuffle(repmat(jit,1, ceil(p.numTrials/numel(jit))));
t.trialJit = trialJit(1:p.numTrials);

t.runDur = t.trialDur*p.numTrials + sum(t.trialJit);
% trialStartTimes = (0:t.trialDur:t.trialDur*p.numTrials-t.trialDur) + cumsum([0 t.trialJit(1:end-1)]);
% 
% targetStartTimes = [];
% cueStartTimes = [];
% 
% % calculate theoretical timing of events
% for n = 1:p.numTrials
%     targetTimes = [0 t.targetSOA];
%     cueTimes = [targetTimes(1)-t.cueTargetSOA targetTimes(2)+t.cueTargetSOA];
%     
%     targetTimes = targetTimes + trialStartTimes(n) + t.trialJit(n);
%     targetStartTimes = [targetStartTimes; targetTimes'];
%     
%     cueTimes = cueTimes + trialStartTimes(n) + t.trialJit(n);
%     cueStartTimes = [cueStartTimes; cueTimes'];
% end
% 
% t.targetStartTimes = targetStartTimes;
% t.cueStartTimes = cueStartTimes;
% t.targetEndTimes = targetStartTimes + t.targetDur;
% p.numTargets = numel(targetStartTimes);

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
surroundGrating = NaN(p.numTrials*2, p.surroundSize, p.surroundSize);

for nTrial = 1:p.numTrials
    center = (sin(p.freq*2*pi/p.centerSize*(Xc.*sin(p.orientation*(pi/180))+Yc.*cos(p.orientation*(pi/180)))-p.phase(nTrial,1)));
    centerGratings(nTrial,:,:) = (center .* centerGaussian);
    
    surround = (sin(p.freqSurround*2*pi/p.surroundSize*(Xs.*sin(p.orientation*(pi/180))+Ys.*cos(p.orientation*(pi/180)))-p.phase(nTrial,2)));
    surroundGrating(nTrial,:,:) = (surround .* surroundGaussian);
end

%% WINDOW SETUP
[window,rect] = Screen('OpenWindow', max(screens), p.grey,[],[],[],[],16);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
% load('MyGammaTable.mat');
% Screen('LoadNormalizedGammaTable', window, repmat(gammaTable, [1 3]));
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
% Draw some text to the screen first outside of the experimental loop:

% Experiment setup
space = zeros(1,256); space(KbName('Space')) = 1;
PsychHID('KbQueueCreate',deviceNumber, space);
PsychHID('KbQueueStart', deviceNumber);


welcomeText = ['Welcome!' '\n' '\n'...
    '' '\n' '\n'...
    'On every trial two center targets are presented next to fixation, on the left and right.' '\n' '\n' ...
    'Your task is to detect a change in contrast in either one of these center targets. ' '\n' '\n' ...
    'An auditory cue (pre-cue) will indicate when a change will occur.' '\n' '\n' ...
    'On most trials the center targets will be accompanied by a surround stimulus, while on some the surround stimulus is absent.' '\n' '\n' ...
    'Additionally, on some trials there will be no auditory cue.' '\n' '\n' ...
    'When asked to, report which target the change in contrast occured.' '\n' '\n' ...
    'Be sure to always maintain steady fixation on the green dot! ' '\n' '\n' '\n' ...
    'Press the UP arrow to continue.' '\n' '\n' ];

DrawFormattedText(window, welcomeText, 'center', 'center', 255);
Screen('Flip', window);
welcomeStart = GetSecs;

% play tones for participant before experiment starts
for nCue=1:size(cueTones,1)
   playSound(pahandle, cueTones(nCue,:)*soundAmp);
   WaitSecs(1);
end

% while 1
%     [pmButton, ~] = PsychPowerMate('Get', powermate);
%     if pmButton == 1
%         break;
%     end
% end

% Make sure we can press esc to quit the experiment
startKey = zeros(1,256); startKey(KbName({'ESCAPE'})) = 1;
PsychHID('KbQueueCreate', deviceNumber, startKey);

% Preallocate some variables
data.estimatedContrast = NaN(p.numTrials, 1);
data.differenceContrast = NaN(p.numTrials, 1);
data.responseTime = NaN(p.numTrials, 1);

% Make surround and center textures before trial loop
surroundStimulus = nan(p.numTrials,1);
centerStimulus = nan(p.numTrials,1);
targetStimulus = nan(p.numTrials,1);

for nTrial = 1:p.numTrials
    surroundTexture(:,:,1) = squeeze(surroundGrating(nTrial,:,:)) * (p.surroundContrast* p.grey ) + p.grey;
%     surroundTexture(:,:,2) = surroundTransparencyMask;
    surroundStimulus(nTrial) = Screen('MakeTexture', window, surroundTexture);
    
    centerTexture(:,:,1) = squeeze(centerGratings(nTrial,:,:)) * ( p.fixedStimContrast * p.grey ) + p.grey;
    centerTexture(:,:,2) = centerTransparencyMask;
    centerStimulus(nTrial) = Screen('MakeTexture', window, centerTexture);

    if mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t1 contrast changes (left target)
        targetTexture(:,:,1) = squeeze(centerGratings(nTrial,:,:)) * ( p.targContrast * p.grey ) + p.grey; 
    elseif mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t2 contrast changes (right target)
        targetTexture(:,:,1) = squeeze(centerGratings(nTrial,:,:)) * ( p.targContrast * p.grey ) + p.grey; 
    end
    targetTexture(:,:,2) = centerTransparencyMask;
    targetStimulus(nTrial) = Screen('MakeTexture', window, targetTexture);   
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
trialTimes = zeros(p.numTrials,6); % [startTrial preCueTime targetTime endTrial]

triggersBase = 1:6;
triggerTimes = nan(p.numTrials,length(triggersBase));
triggers = nan(p.numTrials,length(triggersBase));

for nTrial = 1:p.numTrials
    for nEvent = 1:length(triggersBase)
        triggerChar = [num2str(nTrial) num2str(triggersBase(nEvent))];
        triggers(nTrial,nEvent) = str2double(triggerChar);
    end
end

expStart = GetSecs; % baseline experiment start time
welcomeTime = GetSecs - welcomeStart;
t.welcomeTime = welcomeTime;

for nTrial = 1:p.numTrials
    
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
        shift = randsample(1:p.numTrials,1);
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,4))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,4))
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
        shift = randsample(1:p.numTrials,1);      
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,4))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,4))
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end
    
    % Flickering Target (dur = t.targetDur)
    for nShift = 1:(t.targetDur/t.flicker)
        shift = randsample(1:p.numTrials,1);
        % Draw surroundStimulus if not baseline condition
        if p.trialEvents(nTrial,1) ~= 5 || p.trialEvents(nTrial,1) ~= 6  
            Screen('DrawTexture', window, surroundStimulus(nTrial), [], CenterRectOnPoint([0 0 p.surroundSize p.surroundSize], centerX, centerY), p.trialEvents(nTrial,5))
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)')
            Screen('FillOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)')
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(1), centerY)', p.gapSize, p.gapSize)
            Screen('FrameOval', window, p.grey, CenterRectOnPoint([0 0 p.centerSize+p.gapSize p.centerSize+p.gapSize], patch(2), centerY)', p.gapSize, p.gapSize)
        end
        
        % Determine which target to change
        if mod(p.trialEvents(nTrial,1),2) ~= 0 % if stimConfig is odd, t1 contrast changes (left target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, targetStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,4))
            % Draw centerStimulus2
            Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,4)) 
        elseif mod(p.trialEvents(nTrial,1),2) == 0 % if stimConfig is even, t2 contrast changes (right target)
            % Draw centerStimulus1
            Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,4))
            % Draw centerStimulus2
            Screen('DrawTexture', window, targetStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,4)) 
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
        shift = randsample(1:p.numTrials,1);
        % Draw centerStimulus1
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(1), centerY), p.trialEvents(nTrial,4))
        % Draw centerStimulus2
        Screen('DrawTexture', window, centerStimulus(shift), [], CenterRectOnPoint([0 0 p.centerSize p.centerSize], patch(2), centerY), p.trialEvents(nTrial,4))
        % Draw Fixation
        Screen('FillOval', window, green, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
        Screen('Flip', window);
        WaitSecs(t.flicker);
    end
    
    % Get Response
    
%     % Set up button press
%     PsychHID('KbQueueStart', deviceNumber);
%     PsychHID('KbQueueFlush');
%        
%     % Report: change fixation to allow user to report left/right
%     
    Screen('FillOval', window, blue, [centerX-p.outerFixation centerY-p.outerFixation centerX+p.outerFixation centerY+p.outerFixation])
    Screen('Flip', window);
    GetClicks;
%     
%     startTrial = GetSecs - expStart; % get the start time of each trial
%     
%     %Clicks / scroll
% 
%     % check if esc button has been pressed
%     [keyIsDown, keyCode] = PsychHID('KbQueueCheck', deviceNumber); %check response
%     key = find(keyCode);
%     if key == KbName('ESCAPE') % windows = 'esc', mac = 'ESCAPE' If user presses ESCAPE, exit the program.
%         Screen('LoadNormalizedGammaTable', window, OriginalCLUT);
%         Screen('CloseAll');
%         ListenChar(1); % % Go back to unsuppressed mode
%         FlushEvents('keyDown', deviceNumber);
%         error('User exited program.');
%     end

    % Present center fixation; get ready for next trial
    Screen('Flip', window);
    WaitSecs(t.iti);
 
    %%% Rest period
    if nTrial == p.numTrialsPerBreak*nBreak && nTrial ~= p.numTrialsPerBlock*nBlock
        rest = GetSecs;
        
        restText = 'You can take a short break now, or press the UP arrow to continue.';
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);
        
        nBreak = nBreak+1; 
        
        pmButtonBreak = 0;
         
%         while 1
%             [pmButtonBreak, a] = PsychPowerMate('Get', powermate);
%             if pmButtonBreak == 1
%                 break;
%             end
%         end
        
        t.restTime = (GetSecs-rest)/60;
        GetClicks;
    elseif nTrial == p.numTrialsPerBlock*nBlock
        rest = GetSecs;
        
        restText = ['Block ' num2str(nBlock) ' of ' num2str(p.numBlocks) ' completed! You can take a short break now, ' '' '\n' ...
            'or press the UP arrow to continue' '\n' '\n' ];
        DrawFormattedText(window, restText, 'center', 'center', white);
        Screen('Flip', window);
        
        nBlock = nBlock+1;
        nBreak = nBreak+1; 
        
        pmButtonBreak = 0;
         
%         while 1
%             [pmButtonBreak, a] = PsychPowerMate('Get', powermate);
%             if pmButtonBreak == 1
%                 break;
%             end
%         end
        
        t.restTime = (GetSecs-rest)/60;
        GetClicks;
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

%% timeEvents
t.trialTimes = trialTimes;
timeEvents = nan(p.numTrials, size(t.trialTimes,2)-1); %[cueLeadTime cueTargetSOA targetSOA cueTargetSOA responseTime]

for nEvent = 1:size(t.trialTimes,2)-1
    for nTrial = 1:p.numTrials
        timeEvents(nTrial, nEvent) = t.trialTimes(nTrial,nEvent+1) - t.trialTimes(nTrial,nEvent);
    end
end

t.timeEvents = timeEvents;

% Check timing
% [t.cueStartTimes(1:2:end)+t.welcomeTime t.trialTimes(:,2) t.targetStartTimes(1:2:end)+t.welcomeTime t.trialTimes(:,3) ...
% t.targetStartTimes(2:2:end)+t.welcomeTime t.trialTimes(:,4) t.cueStartTimes(2:2:end)+t.welcomeTime t.trialTimes(:,5)]
%% SAVE OUT THE DATA FILE
cd(dataDir);
theData(p.runNumber).t = t;
theData(p.runNumber).p = p;
theData(p.runNumber).data = data;
eval(['save vTA_surrSuppression_', p.subject, '.mat theData'])

cd(expDir);
