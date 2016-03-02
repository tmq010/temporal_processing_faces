

function [IsQuit, Performance] = test_aaaab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is part of the code to run temporal tuning experiment
% written by Moqian Tian, July 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define parameters

subjID = [];
whichFreq = [];
whichGender = [];
whichOrient = [];
blockDur = [];
totalTargNum = [];
deContrast = [];


while isempty(subjID)
    subjID = input('Please enter your name: ','s');
end

while isempty(whichFreq)
    whichFreq = input('which frequency?: ');
end


while isempty(whichGender)
    whichGender = input('which face gender? 1 female, 2 male: ');
end


while isempty(whichOrient)
    whichOrient = input('which orientation? 1 up, 2 down: ');
end


while isempty(blockDur)
    blockDur = input('How many seconds do you want a block be?: ');
end

while isempty(totalTargNum)
    totalTargNum = input('How many targets do you want in a block?: ');
end

while isempty(deContrast)
    deContrast = input('Decrease the contrast (enter from 0 to 1): ');
end


load('images.mat'); %load the presaved image files

%% here we start
isDone = 0;

% 'endpoint1'
Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'VisualDebugLevel', 1);

% Screen('Preference', 'Verbosity', 1);
% if nargin ~= 1
%     error('Please input all the parameters');
% end
AssertOpenGL; % check for Opengl compatibility, abort otherwise (to verify PTB3)
Screen('Preference', 'SkipSyncTests',1);


%% General Experiment parameters
pixRect = [0 0 206 281]*2;

numBlocks = 1;
sKey = KbName('space');
tarKey = KbName('r');  %% press 'f' when target appears
escKey = KbName('escape');

rootDir = pwd;
fmt = 'bmp';

BeginTime=clock;
Filefix = [subjID,'_freq_',num2str(whichFreq),'_gender_',num2str(whichGender),...
    '_orient_',num2str(whichOrient),'_',num2str(BeginTime(1)), ...
    num2str(BeginTime(2),'%02d'),num2str(BeginTime(3),'%02d'),...
    '_H',int2str(BeginTime(4)),'M',int2str(BeginTime(5)),'S',int2str(BeginTime(6))];


%% Special Experiment parameters
% whichFreq = 2;
% whichGender = 1;
% whichOrient = 1;
%
% blockDur = 12; % total duration of a block in sec
% totalTargNum = 3; % total number of targets in a block

itemsPerBlock = 25; % each gender has 25 examplars
imagesPerBlock = floor(blockDur*whichFreq); % higher freq needs more images
itemDur =  1/whichFreq; % duration of each image
refFrame = 1/60; % single frame duration for 60hz monitor
framesPerItem = floor(itemDur/refFrame); % number of frames per image
% firstFixationDur = 1;
responseWindow = 2; % only count key presses 2s after a target

% Initial screen
screens = Screen('Screens');
screenNumber = max(screens);
screenRect  = Screen(screenNumber, 'rect');

stimuliDir = char([rootDir '/Stimuli']);
dataDir = char([rootDir, '/data/']);

% 'endpoint2'

IsQuit = 0;
try
    
    %     PsychDebugWindowConfiguration(0, [,0.5]);
    
    %% Open a fullscreen window
    backgroundcolor= [138 138 138];
    [window, screenRect] = Screen('OpenWindow',screenNumber,backgroundcolor,screenRect);
%     refresh = Screen('GetFlipInterval', window);
    
%     maskScrPtr = Screen('OpenOffscreenWindow',window,  [0 0 0], pixRect);
    fixation = Screen('OpenOffscreenWindow',window, [138 138 138], screenRect);
    Screen( 'FillOval',fixation, [0 0 0], CenterRect([0 0 12 12], screenRect));
    fixation2 = Screen('OpenOffscreenWindow',window, [138 138 138], screenRect);
    Screen( 'FillOval',fixation2, [0 0 0], CenterRect([0 0 8 8], screenRect));
%     blank = Screen('OpenOffscreenWindow',window,  [138 138 138], screenRect);
    
    
    %% put images in the Screen using 'PutImage'
    cd(stimuliDir);
    for gender=1:2
        cd(['G',num2str(gender)]);
        StimuliFile = dir('*.bmp');
        [numitems, junk] = size(StimuliFile);
        if numitems~=itemsPerBlock, error('Not the right number of items.'); end
        [itemlist{1:numitems}] = deal(StimuliFile.name);
        
        for theitem = 1:itemsPerBlock
            % 		filename = itemlist{theitem};
            % 		[imgArray] = imread(filename, fmt);
            
            for theframe = 1:framesPerItem
                
                % create enough offscreen windows for each picture in the experiment
                offScrPtr(gender,theitem,theframe) = Screen('OpenOffscreenWindow',window,  [138 138 138], pixRect*1.1);
                
                %vary the level of image contrast by a sinusoidal function:
                %0% contrast is gray background, 1% is full image
                a=0:1/(framesPerItem-1):1;
                contrast = sin(a(theframe)*pi); % levels of contrast
                theimage = (imgArray{gender,theitem}-138)*contrast*deContrast+138;
                
                if whichOrient==1
                    % put images into Screen background memory
                    Screen('PutImage',offScrPtr(gender, theitem, theframe),theimage, pixRect*1.1);
                elseif whichOrient==2
                    Screen('PutImage',offScrPtr(gender, theitem, theframe),flipdim(theimage,1), pixRect*1.1);
                end
                
            end
            
        end
        cd(stimuliDir);
    end
    
%     'endpoint5'
    
    % %% read mask image (in case we wanna put masks in between)
    % 	cd(maskDir);
    % 	StimuliFile = dir('*.jpg');
    % 	[numitems junk] = size(StimuliFile);
    % 	[itemlist{1:numitems}] = deal(StimuliFile.name);
    %
    % 	filename = itemlist{1};
    % 	[imgArray] = imread(filename, fmt);
    % 	Screen('PutImage',maskScrPtr, imgArray, pixRect);
    % 	cd(rootDir);
    
    %% define the design matrix
    %% design(numBlocks, imagesPerBlock, 5):
    % 1: image code: randomly choose from 25 exemplars with constraints
    % 2: target code: 1 for target, 0 for non target
    % 3: subject response: 1 for response, 0 for non response
    % 4: image time: absolute time of each image (at the end) measured using
    % function getSecs
    % 5: time when key press: absolute time when a response key is pressed
    % 6: hit
    % 7: false alarm
    % 8: response time: time elapsed between key press and the target
    
    
%     'endpoint6'
    
    
    design = zeros(numBlocks, imagesPerBlock, 8);
    
    
    
    for i = 1:numBlocks
        
        
        % randomly choose target+1 images
        repeatlist = randsample(1:25,1); %define which faces should be repeating
        targlist = randsample(setdiff(1:25,repeatlist),totalTargNum); %define which faces should be the target
        
        %insert target in the images
        allList= round(responseWindow*whichFreq+1):1:round(imagesPerBlock-responseWindow*whichFreq-1);
        targ_tmp = randsample(1:floor(length(allList)/totalTargNum),totalTargNum)*totalTargNum;
        targ = targ_tmp+round(responseWindow*whichFreq);
        targs=sort(targ);
        % always have xx targets within a block
        % don't allow target to appear in the last 2 sec of a block or the
        % first 2 sec of a block. 2 sec should be substituted by how
        % predefined responseWindow
        % targets are never next to each other
        
        design(i,targ,2) = 1; %indexing which one is target
        
        tmplist{1} = [repmat(repeatlist,1,targs(1)-1),targlist(1)];
        for j=2:length(targs)
            tmplist{j} = [repmat(repeatlist,1,targs(j)-targs(j-1)-1),targlist(j)];
        end
            tmplist{j+1} = [repmat(repeatlist,1,imagesPerBlock-targs(end))];
        
        tmplist1=[];
        for k=1:length(targs)+1
        tmplist1 =[tmplist1 tmplist{k}];
        end

        design(i,:,1) =tmplist1;
        
    end
    
    
    
%     'endpoint7'
    
    %% start display
    HideCursor;
    experimentStart = GetSecs;
    %Cue word to remind task and press 'space' to start
    
    for theBlock=1:numBlocks
        
        Screen('TextSize',window,30);
        Screen('TextFont',window,'Arial');
        Screen('TextStyle', window ,1);
        
        Screen('TextColor', window, [0 0 0]);
        Screen('DrawText',window,'Now we will start the test block. Press space to start.' ,...
            (screenRect(3)/2-275),screenRect(4)/2-150);
        Screen('DrawText',window,'Press r when two same faces appear one after another.' ,...
            (screenRect(3)/2-500),screenRect(4)/2-100);
        Screen('Flip',window);
        
        
        while KbCheck, end;
        while 1
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyIsDown && keyCode(sKey)
                break; % start only when start key is pressed
            elseif keyIsDown && keyCode(escKey)
                IsQuit=1;
                break;
            end
        end
        
        %             %% present a fixation screen before stimulus
        %             Screen('CopyWindow', fixation, window, [], screenRect);
        %             startTime = GetSecs;
        %             Screen('Flip',window); % display first fixation
        %
        %             while ~IsQuit & ((GetSecs - startTime)<firstFixationDur)
        %                 [keyIsDown,secs,keyCode] = KbCheck;
        %                 if keyIsDown & keyCode(escKey),
        %                         IsQuit=1;
        %                     break;
        %                 end
        %             end
        % 'endpoint8'
        
        keyWasUp = 1;
        endpointtest(i,1) = GetSecs;
        %% display main experiment
        for i = 1:imagesPerBlock
            
            
            %% present images
            
            thepic = design(theBlock,i,1);
            pixRect_scale = pixRect*randsample([0.9:0.001:1.1],1);
            %scale image randomly between 90% and 110%
            for theframe = 1:framesPerItem
                Screen('CopyWindow', offScrPtr(whichGender,thepic,theframe), window, [],...
                    CenterRectOnPoint(pixRect_scale, screenRect(3)/2,screenRect(4)/2));
                presStart = GetSecs;
                Screen('Flip',window);  % display learning pictures
                
                % check whether this frame should break
                while ~IsQuit && ((GetSecs - presStart)<refFrame)
                    [keyIsDown,secs,keyCode] = KbCheck;
                    
                    if keyIsDown,
                        if keyCode(escKey),
                            IsQuit=1;
                            break;
                        end
                        if keyCode(tarKey) && keyWasUp,
                            design(theBlock,i,3) = 1; % detect a response
                            design(theBlock,i,5) = GetSecs; % time stamp the response
                        end
                        keyWasUp = 0;
                    else
                        keyWasUp = 1;
                    end
                    
                end
                
            end
            
            if IsQuit==1
                break;
            end
            
            design(theBlock,i, 4) = GetSecs; % time stamp of each image
            
        end
        if IsQuit==1
            break;
        end
        endpointtest2 = GetSecs;
    end
    
    % 'endpoint9'
    %% calculate performance
    
    rsp = find(design(1,:,3)==1);
    
    if isempty(rsp)
        data.hit = 0;
        data.fa = 0;
        data.rt = 0;
    elseif length(rsp)>=1
        for i=1:length(rsp)
            for j=0:ceil(responseWindow*whichFreq)
                if rsp(i) < ceil(responseWindow*whichFreq)%if the key press is before the first possible target
                    design(1,rsp(i),7)= 1; %then it's a false alarm
                    design(1,rsp(i),6)= 0;
                    design(1,rsp(i),8)= 0;
                    
                elseif design(1,rsp(i)-j,2)==1 %if key press within 2s window of a target
                        design(1,rsp(i),6) = 1; % then it's a hit
                        design(1,rsp(i),8) = design(1,rsp(i),5)-design(1,rsp(i)-j-1,4);
                        %only calculate RT for correct trials
                        
                elseif imagesPerBlock-rsp(i)<=j %if key press in the last 2s of a block
                        design(1,rsp(i),6) = 0; % then it's a false alarm
                        design(1,rsp(i),8) = 0;
                        design(1,rsp(i),7)= 1;
                end
            end
        end
        
        %     'endpoint10'
        
        % below is to eliminate consecutive key presses between two targets
        starg = sort(targ);
        for k=1:length(starg)-1
            consc=find(design(1,starg(k):1:starg(k+1)-1,6)==1);
            if length(consc)>1
                design(1,starg(k)+consc(2:end)-1,6)=0;
                design(1,starg(k)+consc(2:end)-1,7)=1;
                design(1,starg(k)+consc(2:end)-1,8)=0;
            end
        end

        consc2=find(design(1,starg(end):end,6)==1);
        if length(consc2)>1
            design(1,starg(end)+consc2(2:end)-1,6)=0;
            design(1,starg(end)+consc2(2:end)-1,7)=1;
            design(1,starg(end)+consc2(2:end)-1,8)=0;
        end
        
        %     'endpoint11'
        % false alarm: key presses that are not in 2s window after a target
        design(1,:,7) = design(1,:,3);
        design(1,design(1,:,6)==1,7)=0;
        
        %     'endpoint12'
        % compute summary performance
        % hit and fa are in proportion to total target number
        
        data.hit = length(find(squeeze(design(:,:,6))==1))/totalTargNum;
        data.fa = length(find(squeeze(design(:,:,7))==1))/totalTargNum;
        data.rt = sum(squeeze(design(:,design(:,:,6)==1,8)))/length(rsp);
        
    end
    
    Performance = data;
    %% close up and save data
    Screen('CloseAll');
    cd(rootDir);
    if IsQuit == 1
        disp('ESC is pressed to abort the program.');
        return;
    end
    Screen('CloseAll');
    ShowCursor;
    cd(dataDir);
    save([Filefix,'.mat']);
    cd(rootDir)
    data
catch
    Screen('CloseAll');
    ShowCursor;
    cd(rootDir);
    disp('program error!!!.');
    
    
end



