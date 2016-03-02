%% this piece of code is to generate averaged performance across all blocks

clx

    
    dataFile = dir('*.mat');
	[numfiles, junk] = size(dataFile);
	[itemlist{1:numfiles}] = deal(dataFile.name);
    perfm.hit = zeros(numfiles*32,1);
    perfm.fa = zeros(numfiles*32,1);
    param.gender = zeros(numfiles*32,1);
    param.freq = zeros(numfiles*32,1);
    param.orient = zeros(numfiles*32,1);
    
	for thefile = 1:numfiles
		filename = itemlist{thefile};
		load(filename);
        perfm.hit(1+(thefile-1)*32:thefile*32) = data.hit; %concatenate different blocks into a long vector
        perfm.fa(1+(thefile-1)*32:thefile*32) = data.fa;
        param.gender(1+(thefile-1)*32:thefile*32) = whichGender;
        param.freq(1+(thefile-1)*32:thefile*32) = whichFreq;
        param.orient(1+(thefile-1)*32:thefile*32) = repmat(whichOrient,1,32);
    end
    
    %% average accuracy for each frequency across orientation and gender
    
    acc_avg(1)=(mean(perfm.hit(param.freq==1))-mean(perfm.fa(param.freq==1)))*2;
    acc_avg(2)=(mean(perfm.hit(param.freq==2))-mean(perfm.fa(param.freq==2)))*2;
    acc_avg(3)=(mean(perfm.hit(param.freq==3))-mean(perfm.fa(param.freq==3)))*2;
    acc_avg(4)=(mean(perfm.hit(param.freq==4))-mean(perfm.fa(param.freq==4)))*2;
    acc_avg(5)=(mean(perfm.hit(param.freq==5))-mean(perfm.fa(param.freq==5)))*2;
    acc_avg(6)=(mean(perfm.hit(param.freq==6))-mean(perfm.fa(param.freq==6)))*2;
    acc_avg(7)=(mean(perfm.hit(param.freq==7.5))-mean(perfm.fa(param.freq==7.5)))*2;
    acc_avg(8)=(mean(perfm.hit(param.freq==8.57))-mean(perfm.fa(param.freq==8.57)))*2;
    
    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 0.4])
    plot(acc_avg)
    FreqLabels={'1hz','2hz','3hz','4hz','5hz','6hz','7.5hz','8.57hz'};
    title('averaged accuracy across orientation & gender','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    
    %% average accuracy for each frequency across gender, separte for different orientations
    % hit_up: each row is a frequency, each column is a trial
    for i=1:numfiles*32
        if param.freq(i)==1&&param.orient(i)==1
            hit_up(1,i)=perfm.hit(i);
            fa_up(1,i)=perfm.fa(i);
        elseif param.freq(i)==2&&param.orient(i)==1
            hit_up(2,i)=perfm.hit(i);
            fa_up(2,i)=perfm.fa(i);
        elseif param.freq(i)==3&&param.orient(i)==1
            hit_up(3,i)=perfm.hit(i);
            fa_up(3,i)=perfm.fa(i); 
        elseif param.freq(i)==4&&param.orient(i)==1
            hit_up(4,i)=perfm.hit(i);
            fa_up(4,i)=perfm.fa(i);
        elseif param.freq(i)==5&&param.orient(i)==1
            hit_up(5,i)=perfm.hit(i);
            fa_up(5,i)=perfm.fa(i);
        elseif param.freq(i)==6&&param.orient(i)==1
            hit_up(6,i)=perfm.hit(i);
            fa_up(6,i)=perfm.fa(i);
        elseif param.freq(i)==7.5&&param.orient(i)==1
            hit_up(7,i)=perfm.hit(i);
            fa_up(7,i)=perfm.fa(i);
        elseif param.freq(i)==8.57&&param.orient(i)==1
            hit_up(8,i)=perfm.hit(i);
            fa_up(8,i)=perfm.fa(i);
        elseif param.freq(i)==1&&param.orient(i)==2
            hit_inv(1,i)=perfm.hit(i);
            fa_inv(1,i)=perfm.fa(i);
        elseif param.freq(i)==2&&param.orient(i)==2
            hit_inv(2,i)=perfm.hit(i);
            fa_inv(2,i)=perfm.fa(i);
        elseif param.freq(i)==3&&param.orient(i)==2
            hit_inv(3,i)=perfm.hit(i);
            fa_inv(3,i)=perfm.fa(i); 
        elseif param.freq(i)==4&&param.orient(i)==2
            hit_inv(4,i)=perfm.hit(i);
            fa_inv(4,i)=perfm.fa(i);
        elseif param.freq(i)==5&&param.orient(i)==2
            hit_inv(5,i)=perfm.hit(i);
            fa_inv(5,i)=perfm.fa(i);
        elseif param.freq(i)==6&&param.orient(i)==2
            hit_inv(6,i)=perfm.hit(i);
            fa_inv(6,i)=perfm.fa(i);
        elseif param.freq(i)==7.5&&param.orient(i)==2
            hit_inv(7,i)=perfm.hit(i);
            fa_inv(7,i)=perfm.fa(i);
        elseif param.freq(i)==8.57&&param.orient(i)==2
            hit_inv(8,i)=perfm.hit(i);
            fa_inv(8,i)=perfm.fa(i);
        end
    end
    
    %% average
    for j=1:8
    acc_up(j,1)=(size(find(hit_up(j,:)==1),2)-size(find(fa_up(j,:)==1),2))/numfiles;
    acc_inv(j,1)=(size(find(hit_inv(j,:)==1),2)-size(find(fa_inv(j,:)==1),2))/numfiles;
    end
    acc_diff = acc_up-acc_inv;
    
    
    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 1.2]);
    title('ran');
    subplot(3,1,1);
    plot(acc_up);
    title('averaged accuracy for upright faces','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off');
    subplot(3,1,2)
    plot(acc_inv)
    title('averaged accuracy for inverted faces','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    subplot(3,1,3)
    plot(acc_diff)
    title('averaged accuracy for difference between upright and inverted','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
