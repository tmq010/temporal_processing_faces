

    function [IsQuit, Performance] = test_abcdd_block_try(whichOrient)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% This is part of the code to run temporal tuning experiment
    % written by Moqian Tian, July 2014.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define parameters

    subjID = [];



    while isempty(subjID)
        subjID = input('Please enter your name: ','s');
    end

    trialDur = 5;


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
    refFrame = 1/60; % single frame duration for 60hz monitor

    sKey = KbName('space');
    tarKey = KbName('f');  %% press 'f' when target appears
    nontarKey = KbName('j');
    escKey = KbName('escape');

    rootDir = pwd;
    fmt = 'bmp';

    BeginTime=clock;
    Filefix = [subjID,'_',num2str(BeginTime(1)), ...
        num2str(BeginTime(2),'%02d'),num2str(BeginTime(3),'%02d'),...
        '_H',int2str(BeginTime(4)),'M',int2str(BeginTime(5)),'S',int2str(BeginTime(6))];


    %% Special Experiment parameters

%     numTrials = 64;
%     param_tmp= zeros(3,numTrials);
%     param_tmp(1,:) = repmat([1 2 3 4 5 6 7.5 8.57],1,8);
%     param_tmp(2,:) = [zeros(1,16)+1 ones(1,16)+1 zeros(1,16)+1 ones(1,16)+1];
    numTrials = 32;
    param_tmp= zeros(3,numTrials);
    param_tmp(1,:) = repmat([1 2 3 4 5 6 7.5 8.57],1,4);
    param_tmp(2,:) = [zeros(1,8)+1 ones(1,8)+1 zeros(1,8)+1 ones(1,8)+1];
    param_tmp(3,:) = [zeros(1,16) ones(1,16)];
    param_tmp2=shuffle(param_tmp,2);
    param_tmp3=shuffle(param_tmp2,2);

    whichFreq = param_tmp3(1,:);
    whichGender = param_tmp3(2,:);
    totalTargNum = param_tmp3(3,:);

    itemsPerBlock = 25; % each gender has 25 examplars
    imagesPerTrial = floor(trialDur*whichFreq); % higher freq needs more images
    itemDur =  1./whichFreq; % duration of each image

    framesPerItem = floor(itemDur/refFrame); % number of frames per image



    % Initial screen
    screens = Screen('Screens');
    screenNumber = max(screens);
    screenRect  = Screen(screenNumber, 'rect');

    stimuliDir = char([rootDir '/Stimuli']);
    dataDir = char([rootDir, '/data/']);

    'endpoint2'

    IsQuit = 0;
    try

PsychDebugWindowConfiguration(0, [,0.5]);

        %% Open a fullscreen window
        backgroundcolor= [138 138 138];
        [window, screenRect] = Screen('OpenWindow',screenNumber,backgroundcolor,screenRect);
    %     refresh = Screen('GetFlipInterval', window);

    %     maskScrPtr = Screen('OpenOffscreenWindow',window,  [0 0 0], pixRect);
        fixation = Screen('OpenOffscreenWindow',window, [138 138 138], screenRect);
        Screen( 'FillOval',fixation, [0 0 0], CenterRect([0 0 12 12], screenRect));
    %     blank = Screen('OpenOffscreenWindow',window,  [138 138 138], screenRect);




        'endpoint5'

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
        %% design(numTrials, imagesPerTrial, 5):
        % 1: image code: randomly choose from 25 exemplars with constraints
        % 2: target code: 1 for target, 0 for non target
        % 3: subject response: 1 for response, 0 for non response
        % 4: image time: absolute time of each image (at the end) measured using
        % function getSecs
        % 5: time when key press: absolute time when a response key is pressed
        % 6: hit
        % 7: false alarm
        % 8: response time: time elapsed between key press and the target


%         'endpoint6'





        for i = 1:numTrials

            design{i} = zeros(imagesPerTrial(i), 8);
            % to generate a list of 3750 images that  are randomly chosen
            % from 25 exemplars
            thelist = shuffle(shuffle(repmat(shuffle(1:1:itemsPerBlock),1,150)));

            design{i}(1,1) = thelist(1);
            for j = 2:imagesPerTrial(i)
                design{i}(j,1) = thelist(j);
                if design{i}(j,1) == design{i}(j-1,1)
                    design{i}(j,1) = randsample(setdiff(1:itemsPerBlock,design{i}(j-1,1)),1);
                    %if two consecutive images are the same, replace it with a
                    %different one
                end
            end

            if totalTargNum(i)==1
                %insert target in the images
                allList= ceil(imagesPerTrial(i)*0.2+1):1:round(imagesPerTrial(i)*0.8-1);
                targ_tmp = randsample(allList,1);
                targ = targ_tmp;
                % always have xx targets within a block
                % target only appear between 20%-80% trials


                for j=1:imagesPerTrial(i)
                    for k = 1:length(targ)
                        if j==targ(k)
                            design{i}(j,1) = design{i}(j-1,1); % make the target same as the previous image
                        end
                    end
                end


                design{i}(targ,2) = 1; %indexing which one is target

                % if a target is the same as the next image, exchange the next
                % image
                for thetar = 1:length(targ)
                    if design{i}(targ(thetar),1) == design{i}(targ(thetar)+1,1)
                        design{i}(targ(thetar)+1,1) = randsample(setdiff(1:itemsPerBlock,...
                            [design{i}(targ(thetar),1) design{i}(targ(thetar)+2,1)]),1);
                    end
                end     
            end
        end



%         'endpoint7'

        %% start display
        HideCursor;
        experimentStart = GetSecs;

        %Cue word to remind task and press 'space' to start
            Screen('TextSize',window,30);
            Screen('TextFont',window,'Arial');
            Screen('TextStyle', window ,1);

            Screen('TextColor', window, [0 0 0]);
            Screen('DrawText',window,'Now we will start the test block.' ,...
                (screenRect(3)/2-275),screenRect(4)/2-150);
            Screen('DrawText',window,'after each trial, press f to report a target, press j if there is no target' ,...
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
            
        data.hit = zeros(numTrials,1);
        data.fa = zeros(numTrials,1);
        
        %% start experiment 
        
        for theTrial=1:numTrials
        %% putting images for each trial into memory buffer
%     'endpoint7.1'
        gender = whichGender(theTrial);
%             cd(['G',num2str(gender)]);
%             StimuliFile = dir('*.bmp');
%             [numitems, junk] = size(StimuliFile);
%             if numitems~=itemsPerBlock, error('Not the right number of items.'); end
%             [itemlist{1:numitems}] = deal(StimuliFile.name);

            for theitem = 1:itemsPerBlock
                % 		filename = itemlist{theitem};
                % 		[imgArray] = imread(filename, fmt);

                for theframe = 1:framesPerItem(theTrial)

                    % create enough offscreen windows for each picture in the experiment
                    %offScrPtr(theTrial,gender,theitem,theframe) = Screen('OpenOffscreenWindow',window,  [138 138 138], pixRect*1.1);
                    

                    %vary the level of image contrast by a sinusoidal function:
                    %0% contrast is gray background, 1% is full image
                    a=0:1/(framesPerItem(theTrial)-1):1;
                    contrast = sin(a(theframe)*pi); % levels of contrast
                    theimage = (imgArray{gender,theitem}-138)*contrast+138;

                    if whichOrient==1
                        % put images into offscreen buffer
                        offScrPtr(gender,theitem,theframe) = Screen('MakeTexture', window,  theimage);
                        %Screen('DrawTexture',offScrPtr(theTrial,gender, theitem, theframe),theimage, pixRect*1.1);
                    elseif whichOrient==2
                        %Screen('DrawTexture',offScrPtr(theTrial,gender, theitem, theframe),flipdim(theimage,1), pixRect*1.1);
                        offScrPtr(gender,theitem,theframe) = Screen('MakeTexture', window,  flipdim(theimage,1));
                    end

                end
            end
            endpointtest(theTrial,1) = GetSecs;
            
            %% start presenting images
            for i = 1:imagesPerTrial(theTrial)
                thepic = design{theTrial}(i,1);
                pixRect_scale = pixRect*randsample([0.9:0.001:1.1],1);
                %scale image randomly between 90% and 110%
                for theframe = 1:framesPerItem(theTrial)
                    %Screen('DrawTexture', offScrPtr(theTrial,whichGender(theTrial),thepic,theframe), window, [],...
                        %CenterRectOnPoint(pixRect_scale, screenRect(3)/2,screenRect(4)/2));
                    Screen('DrawTexture', window, offScrPtr(whichGender(theTrial),thepic,theframe),[], CenterRectOnPoint(pixRect_scale, screenRect(3)/2,screenRect(4)/2));
                    presStart = GetSecs;
                    Screen('Flip',window);  % display learning pictures
%     'endpoint7.1'
                    % check whether this frame should break
                    while ~IsQuit & ((GetSecs - presStart)<refFrame)
                    [keyIsDown,secs,keyCode] = KbCheck;
                        if keyIsDown & keyCode(escKey),
                            IsQuit=1;
                            break;
                        end
                    end

                end

                if IsQuit==1
                    break;
                end

                design{theTrial}(i, 4) = GetSecs; % time stamp of each image

            end
            if IsQuit==1
                break;
            end
            endpointtest2 = GetSecs;

                %% fixation screen waiting for response
                Screen('CopyWindow', fixation, window, [], CenterRect(pixRect, screenRect));
                Screen('Flip',window); % display fixation
                Response = 0;
                %% record response while fixation
                while ~Response && ~IsQuit  
                    [keyIsDown,secs,keyCode] = KbCheck;
                    if (keyIsDown && ( keyCode(tarKey) || keyCode(nontarKey)||keyCode(escKey)) && length(find(keyCode))==1)
                        Response = 1;
                        switch (find(keyCode))
                            case tarKey
                                design{theTrial}(1,3) = 1;
                            case nontarKey
                                design{theTrial}(1,3) = 0;
                            case escKey
                                IsQuit=1;
                                break;
                        end;
                        while KbCheck;
                        end;
                    end
                end
                
                % calculate performance for each trial
                if size(find(design{theTrial}(:,3)==1),1)==1 && size(find(design{theTrial}(:,2)==1),1)==1
                data.hit(theTrial) = 1;
                elseif size(find(design{theTrial}(:,3)==1),1)==1 && size(find(design{theTrial}(:,2)==1),1)==0
                data.fa(theTrial) = 1;
                end

        end
%         'endpoint9'

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

    catch
        Screen('CloseAll');
        ShowCursor;
        cd(rootDir);
        disp('program error!!!.');


    end



