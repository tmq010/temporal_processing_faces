%% this code is for analyzing the exact timing of stimuli presentation
%%  for the 4s/trial design.

% clx
% 
%     
%     dataFile = dir('*.mat');
% 	[numfiles junk] = size(dataFile);
% 	[itemlist{1:numfiles}] = deal(dataFile.name);
%     
%     
% 	for thefile = 1:numfiles
% 		filename = itemlist{thefile};
% 		load(filename);
%         for i=1:32
%         timeStamps{(thefile-1)*32+i} = design{i}(:,4);
%         end
%         param.freq(1+(thefile-1)*32:thefile*32) = whichFreq;
%     end
numfiles=1;
thefile=numfiles;
        for i=1:32
        timeStamps{(thefile-1)*32+i} = design{i}(:,4);
        end
param.freq(1+(thefile-1)*32:thefile*32) = whichFreq;

        
    ts_1hz = zeros(4,numfiles*32);
    ts_2hz = zeros(8,numfiles*32);
    ts_3hz = zeros(12,numfiles*32);
    ts_4hz = zeros(16,numfiles*32);
    ts_5hz = zeros(20,numfiles*32);
    ts_6hz = zeros(24,numfiles*32);
    ts_7hz = zeros(30,numfiles*32);
    ts_8hz = zeros(34,numfiles*32);
    
    for i = 1:numfiles*32
        if size(timeStamps{i},1)==4
            ts_1hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==8
            ts_2hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==12
            ts_3hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==16
            ts_4hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==20
            ts_5hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==24
            ts_6hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)== 30
            ts_7hz(:,i) = timeStamps{i}(:);
        elseif size(timeStamps{i},1)==34
            ts_8hz(:,i) = timeStamps{i}(:);
        end
    end
    
    for i=1:3
    dur_1hz(i,:) = ts_1hz(i+1,:)-ts_1hz(i,:);
    end
    for i=1:7
    dur_2hz(i,:) = ts_2hz(i+1,:)-ts_2hz(i,:);
    end
    for i=1:11
    dur_3hz(i,:) = ts_3hz(i+1,:)-ts_3hz(i,:);
    end
    for i=1:15
    dur_4hz(i,:) = ts_4hz(i+1,:)-ts_4hz(i,:);
    end
    for i=1:19
    dur_5hz(i,:) = ts_5hz(i+1,:)-ts_5hz(i,:);
    end
    for i=1:23
    dur_6hz(i,:) = ts_6hz(i+1,:)-ts_6hz(i,:);
    end
    for i=1:29
    dur_7hz(i,:) = ts_7hz(i+1,:)-ts_7hz(i,:);
    end
    for i=1:33
    dur_8hz(i,:) = ts_8hz(i+1,:)-ts_8hz(i,:);
    end
    

    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 1.2]);
    subplot(2,4,1)
    boxplot(dur_1hz(dur_1hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('1hz(1000ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([1000 1005]); refline(0,1000)
    
    subplot(2,4,2)
    boxplot(dur_2hz(dur_2hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('2hz(500ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([500 503]); refline(0,500)
    
    subplot(2,4,3)
    boxplot(dur_3hz(dur_3hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('3hz(333.3ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([333 336]);refline(0,1000/3)
    
    subplot(2,4,4)
    boxplot(dur_4hz(dur_4hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('4hz(250ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([250 253]); refline(0,1000/4)
    
    subplot(2,4,5)
    boxplot(dur_5hz(dur_5hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('5hz(200ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([200 205]); refline(0,1000/5)
    
    subplot(2,4,6)
    boxplot(dur_6hz(dur_6hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('6hz(166.7ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([166 170]); refline(0,1000/6)
    
    subplot(2,4,7)
    boxplot(dur_7hz(dur_7hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('7.5hz(133.3ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([132 136]); refline(0,1000/7.5)
    
    subplot(2,4,8)
    boxplot(dur_8hz(dur_8hz~=0)*1000, 'plotstyle','traditional','boxstyle','outline')
    ylabel('8.57hz(116.7ms/image)','Fontsize',16,'FontWeight', 'BOLD');
    set(gca,'XtickLabel','actual ms/image', 'Fontsize',10,'box','off')
    ylim([115 118]);refline(0,1000/8.57)

    