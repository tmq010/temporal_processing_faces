%% this piece of code is to generate averaged performance across blocks
%% the difference between analysis_2.m and analysis.m is that analysis_2.m generate error bars

clx

    
    dataFile = dir('*.mat');
	[numfiles, junk] = size(dataFile);
	[itemlist{1:numfiles}] = deal(dataFile.name);
    perfm.hit = zeros(32,numfiles);
    perfm.fa = zeros(32,numfiles);
    param.gender = zeros(32,numfiles);
    param.freq = zeros(32,numfiles);
    param.orient = zeros(32,numfiles);
    param.targ = zeros(32,numfiles);
    
	for thefile = 1:numfiles
		filename = itemlist{thefile};
		load(filename);
        perfm.hit(:,thefile) = data.hit; 
        perfm.fa(:,thefile) = data.fa;
        param.gender(:,thefile) = whichGender;
        param.freq(:,thefile) = whichFreq;
        param.orient(:,thefile) = repmat(whichOrient,32,1);
        
        for j=1:32
            if isempty(find(design{j}(:,2)==1, 1))==1
               param.targ(j,thefile)=0; % meaning no target in this trial
            elseif isempty(find(design{j}(:,2)==1, 1))==0 &&  find(design{j}(:,2)==1)<length(design{j}(:,2))/2
               param.targ(j,thefile)=1; % meaning an early target
            elseif isempty(find(design{j}(:,2)==1, 1))==0 &&  find(design{j}(:,2)==1)>=length(design{j}(:,2))/2
               param.targ(j,thefile)=2; % meaning a late target 
            end
        end
    end
    
    %% average accuracy for each frequency across orientation and gender
    for thefile = 1:numfiles
    acc_avg(1,thefile)=(mean(perfm.hit(param.freq(:,thefile)==1,thefile))-mean(perfm.fa(param.freq(:,thefile)==1,thefile)))*2;
    acc_avg(2,thefile)=(mean(perfm.hit(param.freq(:,thefile)==2,thefile))-mean(perfm.fa(param.freq(:,thefile)==2,thefile)))*2;
    acc_avg(3,thefile)=(mean(perfm.hit(param.freq(:,thefile)==3,thefile))-mean(perfm.fa(param.freq(:,thefile)==3,thefile)))*2;
    acc_avg(4,thefile)=(mean(perfm.hit(param.freq(:,thefile)==4,thefile))-mean(perfm.fa(param.freq(:,thefile)==4,thefile)))*2;
    acc_avg(5,thefile)=(mean(perfm.hit(param.freq(:,thefile)==5,thefile))-mean(perfm.fa(param.freq(:,thefile)==5,thefile)))*2;
    acc_avg(6,thefile)=(mean(perfm.hit(param.freq(:,thefile)==6,thefile))-mean(perfm.fa(param.freq(:,thefile)==6,thefile)))*2;
    acc_avg(7,thefile)=(mean(perfm.hit(param.freq(:,thefile)==7.5,thefile))-mean(perfm.fa(param.freq(:,thefile)==7.5,thefile)))*2;
    acc_avg(8,thefile)=(mean(perfm.hit(param.freq(:,thefile)==8.57,thefile))-mean(perfm.fa(param.freq(:,thefile)==8.57,thefile)))*2;
    
    hit_avg(1,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==1,thefile))*2;
    hit_avg(2,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==2,thefile))*2;
    hit_avg(3,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==3,thefile))*2;
    hit_avg(4,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==4,thefile))*2;
    hit_avg(5,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==5,thefile))*2;
    hit_avg(6,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==6,thefile))*2;
    hit_avg(7,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==7.5,thefile))*2;
    hit_avg(8,thefile) = nanmean(perfm.hit(param.freq(:,thefile)==8.57,thefile))*2;
    
    fa_avg(1,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==1,thefile))*2;
    fa_avg(2,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==2,thefile))*2;
    fa_avg(3,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==3,thefile))*2;
    fa_avg(4,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==4,thefile))*2;
    fa_avg(5,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==5,thefile))*2;
    fa_avg(6,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==6,thefile))*2;
    fa_avg(7,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==7.5,thefile))*2;
    fa_avg(8,thefile) = nanmean(perfm.fa(param.freq(:,thefile)==8.57,thefile))*2;
    
    end
    
    for freq = 1:8
    std_avg(freq,1) = std(acc_avg(freq,:));
    std_hit(freq,1) = std(hit_avg(freq,:));
    std_fa(freq,1) = std(fa_avg(freq,:));
    end
    
    cilow = mean(acc_avg,2)-1.96*std_avg/sqrt(numfiles);
    cihigh = mean(acc_avg,2)+1.96*std_avg/sqrt(numfiles);
    
    cilow_hit = mean(hit_avg,2)-1.96*std_hit/sqrt(numfiles);
    cihigh_hit = mean(hit_avg,2)+1.96*std_hit/sqrt(numfiles);
    
    cilow_fa = mean(fa_avg,2)-1.96*std_fa/sqrt(numfiles);
    cihigh_fa = mean(fa_avg,2)+1.96*std_fa/sqrt(numfiles);
    
    %% plot average accuracy 
    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 0.4])
    ciplot(cilow, cihigh, [1:8],'r',.5)
    hold on
    ciplot(cilow_hit, cihigh_hit, [1:8],'g',.2)
    ciplot(cilow_fa, cihigh_fa, [1:8],'b',.2)

    
    x=[1:8]';
    h=plot(x,mean(acc_avg,2),'r',x,mean(hit_avg,2),'g',x,mean(fa_avg,2),'b');
    legend('acc','hit','fa');
    for i=1:3
    set(h(i),'LineWidth',1.5);
    end
    
    FreqLabels={'1hz','2hz','3hz','4hz','5hz','6hz','7.5hz','8.57hz'};
%     FreqLabels={'0.9hz','1.9hz','2.8hz','3.7hz','4.7hz','5.7hz','7.1hz','8.1hz'}; % first two blocks of moqian and ran

    title('averaged accuracy across orientation & gender','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    hold off
    
    %% average accuracy for each frequency across gender, separte for different orientations
    % hit_up: each row is a frequency, each column is a trial
        earlyhit_up_1hz = zeros(32,numfiles);
        latehit_up_1hz = zeros(32,numfiles);
        earlyhit_up_2hz = zeros(32,numfiles);
        latehit_up_2hz = zeros(32,numfiles);
        earlyhit_up_3hz = zeros(32,numfiles);
        latehit_up_3hz = zeros(32,numfiles);
        earlyhit_up_4hz = zeros(32,numfiles);
        latehit_up_4hz = zeros(32,numfiles);
        earlyhit_up_5hz = zeros(32,numfiles);
        latehit_up_5hz = zeros(32,numfiles);
        earlyhit_up_6hz = zeros(32,numfiles);
        latehit_up_6hz = zeros(32,numfiles);
        earlyhit_up_7hz = zeros(32,numfiles);
        latehit_up_7hz = zeros(32,numfiles);
        earlyhit_up_8hz = zeros(32,numfiles);
        latehit_up_8hz = zeros(32,numfiles);
        earlyhit_inv_1hz = zeros(32,numfiles);
        latehit_inv_1hz = zeros(32,numfiles);
        earlyhit_inv_2hz = zeros(32,numfiles);
        latehit_inv_2hz = zeros(32,numfiles);
        earlyhit_inv_3hz = zeros(32,numfiles);
        latehit_inv_3hz = zeros(32,numfiles);
        earlyhit_inv_4hz = zeros(32,numfiles);
        latehit_inv_4hz = zeros(32,numfiles);
        earlyhit_inv_5hz = zeros(32,numfiles);
        latehit_inv_5hz = zeros(32,numfiles);
        earlyhit_inv_6hz = zeros(32,numfiles);
        latehit_inv_6hz = zeros(32,numfiles);
        earlyhit_inv_7hz = zeros(32,numfiles);
        latehit_inv_7hz = zeros(32,numfiles);
        earlyhit_inv_8hz = zeros(32,numfiles);
        latehit_inv_8hz = zeros(32,numfiles);
    
    for i=1:numfiles
        if param.orient(:,i)==1
        hit_up_1hz(:,i) = perfm.hit(param.freq(:,i)==1, i);
        fa_up_1hz(:,i) = perfm.fa(param.freq(:,i)==1, i);
        hit_up_2hz(:,i) = perfm.hit(param.freq(:,i)==2, i);
        fa_up_2hz(:,i) = perfm.fa(param.freq(:,i)==2, i);
        hit_up_3hz(:,i) = perfm.hit(param.freq(:,i)==3, i);
        fa_up_3hz(:,i) = perfm.fa(param.freq(:,i)==3, i);
        hit_up_4hz(:,i) = perfm.hit(param.freq(:,i)==4, i);
        fa_up_4hz(:,i) = perfm.fa(param.freq(:,i)==4, i);
        hit_up_5hz(:,i) = perfm.hit(param.freq(:,i)==5, i);
        fa_up_5hz(:,i) = perfm.fa(param.freq(:,i)==5, i);
        hit_up_6hz(:,i) = perfm.hit(param.freq(:,i)==6, i);
        fa_up_6hz(:,i) = perfm.fa(param.freq(:,i)==6, i);
        hit_up_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5, i);
        fa_up_7hz(:,i) = perfm.fa(param.freq(:,i)==7.5, i);
        hit_up_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57, i);
        fa_up_8hz(:,i) = perfm.fa(param.freq(:,i)==8.57, i);

        
        earlyhit_up_1hz(:,i) = perfm.hit(param.freq(:,i)==1&param.targ(:,i)==1, i);
        latehit_up_1hz(:,i) = perfm.hit(param.freq(:,i)==1&param.targ(:,i)==2, i);
        earlyhit_up_2hz(:,i) = perfm.hit(param.freq(:,i)==2&param.targ(:,i)==1, i);
        latehit_up_2hz(:,i) = perfm.hit(param.freq(:,i)==2&param.targ(:,i)==2, i);
        earlyhit_up_3hz(:,i) = perfm.hit(param.freq(:,i)==3&param.targ(:,i)==1, i);
        latehit_up_3hz(:,i) = perfm.hit(param.freq(:,i)==3&param.targ(:,i)==2, i);
        earlyhit_up_4hz(:,i) = perfm.hit(param.freq(:,i)==4&param.targ(:,i)==1, i);
        latehit_up_4hz(:,i) = perfm.hit(param.freq(:,i)==4&param.targ(:,i)==2, i);
        earlyhit_up_5hz(:,i) = perfm.hit(param.freq(:,i)==5&param.targ(:,i)==1, i);
        latehit_up_5hz(:,i) = perfm.hit(param.freq(:,i)==5&param.targ(:,i)==2, i);
        earlyhit_up_6hz(:,i) = perfm.hit(param.freq(:,i)==6&param.targ(:,i)==1, i);
        latehit_up_6hz(:,i) = perfm.hit(param.freq(:,i)==6&param.targ(:,i)==2, i);
        earlyhit_up_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5&param.targ(:,i)==1, i);
        latehit_up_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5&param.targ(:,i)==2, i);
        earlyhit_up_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57&param.targ(:,i)==1, i);
        latehit_up_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57&param.targ(:,i)==2, i);
        
        if param.orient(:,i)==2
        hit_inv_1hz(:,i) = perfm.hit(param.freq(:,i)==1, i);
        fa_inv_1hz(:,i) = perfm.fa(param.freq(:,i)==1, i);
        hit_inv_2hz(:,i) = perfm.hit(param.freq(:,i)==2, i);
        fa_inv_2hz(:,i) = perfm.fa(param.freq(:,i)==2, i);
        hit_inv_3hz(:,i) = perfm.hit(param.freq(:,i)==3, i);
        fa_inv_3hz(:,i) = perfm.fa(param.freq(:,i)==3, i);
        hit_inv_4hz(:,i) = perfm.hit(param.freq(:,i)==4, i);
        fa_inv_4hz(:,i) = perfm.fa(param.freq(:,i)==4, i);
        hit_inv_5hz(:,i) = perfm.hit(param.freq(:,i)==5, i);
        fa_inv_5hz(:,i) = perfm.fa(param.freq(:,i)==5, i);
        hit_inv_6hz(:,i) = perfm.hit(param.freq(:,i)==6, i);
        fa_inv_6hz(:,i) = perfm.fa(param.freq(:,i)==6, i);
        hit_inv_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5, i);
        fa_inv_7hz(:,i) = perfm.fa(param.freq(:,i)==7.5, i);
        hit_inv_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57, i);
        fa_inv_8hz(:,i) = perfm.fa(param.freq(:,i)==8.57, i);
        
        earlyhit_inv_1hz(:,i) = perfm.hit(param.freq(:,i)==1&param.targ(:,i)==1, i);
        latehit_inv_1hz(:,i) = perfm.hit(param.freq(:,i)==1&param.targ(:,i)==2, i);
        earlyhit_inv_2hz(:,i) = perfm.hit(param.freq(:,i)==2&param.targ(:,i)==1, i);
        latehit_inv_2hz(:,i) = perfm.hit(param.freq(:,i)==2&param.targ(:,i)==2, i);
        earlyhit_inv_3hz(:,i) = perfm.hit(param.freq(:,i)==3&param.targ(:,i)==1, i);
        latehit_inv_3hz(:,i) = perfm.hit(param.freq(:,i)==3&param.targ(:,i)==2, i);
        earlyhit_inv_4hz(:,i) = perfm.hit(param.freq(:,i)==4&param.targ(:,i)==1, i);
        latehit_inv_4hz(:,i) = perfm.hit(param.freq(:,i)==4&param.targ(:,i)==2, i);
        earlyhit_inv_5hz(:,i) = perfm.hit(param.freq(:,i)==5&param.targ(:,i)==1, i);
        latehit_inv_5hz(:,i) = perfm.hit(param.freq(:,i)==5&param.targ(:,i)==2, i);
        earlyhit_inv_6hz(:,i) = perfm.hit(param.freq(:,i)==6&param.targ(:,i)==1, i);
        latehit_inv_6hz(:,i) = perfm.hit(param.freq(:,i)==6&param.targ(:,i)==2, i);
        earlyhit_inv_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5&param.targ(:,i)==1, i);
        latehit_inv_7hz(:,i) = perfm.hit(param.freq(:,i)==7.5&param.targ(:,i)==2, i);
        earlyhit_inv_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57&param.targ(:,i)==1, i);
        latehit_inv_8hz(:,i) = perfm.hit(param.freq(:,i)==8.57&param.targ(:,i)==2, i);
        end
    end
    
    
    %% calculation
    
    
    a1=mean(hit_up_1hz); avg_hit_up(1,:) = a1(param.orient(1,:)==1); % exclude blocks that have inverted faces
    a2=mean(hit_up_2hz); avg_hit_up(2,:) = a2(param.orient(1,:)==1); % each row is a frequency
    a3=mean(hit_up_3hz); avg_hit_up(3,:) = a3(param.orient(1,:)==1); % each column is a block
    a4=mean(hit_up_4hz); avg_hit_up(4,:) = a4(param.orient(1,:)==1);
    a5=mean(hit_up_5hz); avg_hit_up(5,:) = a5(param.orient(1,:)==1);
    a6=mean(hit_up_6hz); avg_hit_up(6,:) = a6(param.orient(1,:)==1);
    a7=mean(hit_up_7hz); avg_hit_up(7,:) = a7(param.orient(1,:)==1);
    a8=mean(hit_up_8hz); avg_hit_up(8,:) = a8(param.orient(1,:)==1);
    
    
    b1=mean(hit_inv_1hz); avg_hit_inv(1,:) = b1(param.orient(1,:)==2); % exclude blocks that have upright faces
    b2=mean(hit_inv_2hz); avg_hit_inv(2,:) = b2(param.orient(1,:)==2); % each row is a frequency
    b3=mean(hit_inv_3hz); avg_hit_inv(3,:) = b3(param.orient(1,:)==2); % each column is a block
    b4=mean(hit_inv_4hz); avg_hit_inv(4,:) = b4(param.orient(1,:)==2);
    b5=mean(hit_inv_5hz); avg_hit_inv(5,:) = b5(param.orient(1,:)==2);
    b6=mean(hit_inv_6hz); avg_hit_inv(6,:) = b6(param.orient(1,:)==2);
    b7=mean(hit_inv_7hz); avg_hit_inv(7,:) = b7(param.orient(1,:)==2);
    b8=mean(hit_inv_8hz); avg_hit_inv(8,:) = b8(param.orient(1,:)==2);
    
    c1=mean(fa_up_1hz); avg_fa_up(1,:) = c1(param.orient(1,:)==1); % exclude blocks that have inverted faces
    c2=mean(fa_up_2hz); avg_fa_up(2,:) = c2(param.orient(1,:)==1); % each row is a frequency
    c3=mean(fa_up_3hz); avg_fa_up(3,:) = c3(param.orient(1,:)==1); % each column is a block
    c4=mean(fa_up_4hz); avg_fa_up(4,:) = c4(param.orient(1,:)==1);
    c5=mean(fa_up_5hz); avg_fa_up(5,:) = c5(param.orient(1,:)==1);
    c6=mean(fa_up_6hz); avg_fa_up(6,:) = c6(param.orient(1,:)==1);
    c7=mean(fa_up_7hz); avg_fa_up(7,:) = c7(param.orient(1,:)==1);
    c8=mean(fa_up_8hz); avg_fa_up(8,:) = c8(param.orient(1,:)==1);
    
    
    d1=mean(fa_inv_1hz); avg_fa_inv(1,:) = d1(param.orient(1,:)==2); % exclude blocks that have upright faces
    d2=mean(fa_inv_2hz); avg_fa_inv(2,:) = d2(param.orient(1,:)==2); % each row is a frequency
    d3=mean(fa_inv_3hz); avg_fa_inv(3,:) = d3(param.orient(1,:)==2); % each column is a block
    d4=mean(fa_inv_4hz); avg_fa_inv(4,:) = d4(param.orient(1,:)==2);
    d5=mean(fa_inv_5hz); avg_fa_inv(5,:) = d5(param.orient(1,:)==2);
    d6=mean(fa_inv_6hz); avg_fa_inv(6,:) = d6(param.orient(1,:)==2);
    d7=mean(fa_inv_7hz); avg_fa_inv(7,:) = d7(param.orient(1,:)==2);
    d8=mean(fa_inv_8hz); avg_fa_inv(8,:) = d8(param.orient(1,:)==2);
    
    avg_acc_up = (avg_hit_up - avg_fa_up)*2;
    avg_acc_inv = (avg_hit_inv - avg_fa_inv)*2;
    
    for freq = 1:8
    std_avg_up(freq,1) = std(avg_acc_up(freq,:));
    std_avg_inv(freq,1) = std(avg_acc_inv(freq,:));
    std_hit_up(freq,1) = std(avg_hit_up(freq,:));
    std_hit_inv(freq,1) = std(avg_hit_inv(freq,:));
    std_fa_up(freq,1) = std(avg_fa_up(freq,:));
    std_fa_inv(freq,1) = std(avg_fa_inv(freq,:));
    end
    
    %% acc
    cilow_up = mean(avg_acc_up,2)-1.96*std_avg_up/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_up = mean(avg_acc_up,2)+1.96*std_avg_up/sqrt(size(find(param.orient(1,:)==1),2));
    
    cilow_inv = mean(avg_acc_inv,2)-1.96*std_avg_inv/sqrt(size(find(param.orient(1,:)==2),2));
    cihigh_inv = mean(avg_acc_inv,2)+1.96*std_avg_inv/sqrt(size(find(param.orient(1,:)==2),2));
    
    avg_acc_diff = mean(avg_acc_up,2) - mean(avg_acc_inv,2);
    std_acc_diff = sqrt(std_avg_up.^2/size(find(param.orient(1,:)==1),2)+std_avg_inv.^2/size(find(param.orient(1,:)==2),2));
    
    cilow_diff = avg_acc_diff-1.96*std_acc_diff/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_diff = avg_acc_diff+1.96*std_acc_diff/sqrt(size(find(param.orient(1,:)==1),2));
    
    %% hit
    cilow_hit_up = mean(avg_hit_up,2)-1.96*std_hit_up/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_hit_up = mean(avg_hit_up,2)+1.96*std_hit_up/sqrt(size(find(param.orient(1,:)==1),2));
    
    cilow_hit_inv = mean(avg_hit_inv,2)-1.96*std_hit_inv/sqrt(size(find(param.orient(1,:)==2),2));
    cihigh_hit_inv = mean(avg_hit_inv,2)+1.96*std_hit_inv/sqrt(size(find(param.orient(1,:)==2),2));
    
    avg_hit_diff = mean(avg_hit_up,2) - mean(avg_hit_inv,2);
    std_hit_diff = sqrt(std_hit_up.^2/size(find(param.orient(1,:)==1),2)+std_hit_inv.^2/size(find(param.orient(1,:)==2),2));
    
    cilow_hit_diff = avg_hit_diff-1.96*std_hit_diff/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_hit_diff = avg_hit_diff+1.96*std_hit_diff/sqrt(size(find(param.orient(1,:)==1),2));
    
    %% fa
    cilow_fa_up = mean(avg_fa_up,2)-1.96*std_fa_up/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_fa_up = mean(avg_fa_up,2)+1.96*std_fa_up/sqrt(size(find(param.orient(1,:)==1),2));
    
    cilow_fa_inv = mean(avg_fa_inv,2)-1.96*std_fa_inv/sqrt(size(find(param.orient(1,:)==2),2));
    cihigh_fa_inv = mean(avg_fa_inv,2)+1.96*std_fa_inv/sqrt(size(find(param.orient(1,:)==2),2));
    
    avg_fa_diff = mean(avg_fa_up,2) - mean(avg_fa_inv,2);
    std_fa_diff = sqrt(std_fa_up.^2/size(find(param.orient(1,:)==1),2)+std_fa_inv.^2/size(find(param.orient(1,:)==2),2));
    
    cilow_fa_diff = avg_fa_diff-1.96*std_fa_diff/sqrt(size(find(param.orient(1,:)==1),2));
    cihigh_fa_diff = avg_fa_diff+1.96*std_fa_diff/sqrt(size(find(param.orient(1,:)==1),2));
    
    %% draw graph with the changes only in accuracy
    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 1.2]);
    
    subplot(3,1,1);
    ciplot(cilow_up, cihigh_up, [1:8],'r',.3)
    hold on  
    plot(mean(avg_acc_up,2),'linewidth',1.5,'color','b')
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off');
    ylim([-0.2 1]);
    title('averaged accuracy for upright faces','Fontsize',18,'FontWeight', 'BOLD');
    hold off
    
    subplot(3,1,2);
    ciplot(cilow_inv, cihigh_inv, [1:8],'r',.3)
    hold on  
    plot(mean(avg_acc_inv,2),'linewidth',1.5,'color','b')
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off');
    ylim([-0.2 1]);
    title('averaged accuracy for inverted faces','Fontsize',18,'FontWeight', 'BOLD');
    hold off
    
    subplot(3,1,3);
    ciplot(cilow_diff, cihigh_diff, [1:8],'r',.3)
    hold on  
    plot(mean(avg_acc_diff,2),'linewidth',1.5,'color','b')
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off');
    ylim([-0.5 1]);
    title('inversion effect across frequencies','Fontsize',18,'FontWeight', 'BOLD');
    hold off
    
    %% plot accuracy along with hit and false alarm

    figure('Color',[ 1 1 1],  'units','norm', 'position', [ .1 .1 .4 1.2]);
    
    subplot(3,1,1);
    ciplot(cilow_up, cihigh_up, [1:8],'r',.5)
    hold on
    ciplot(cilow_hit_up, cihigh_hit_up, [1:8],'g',.2)
    ciplot(cilow_fa_up, cihigh_fa_up, [1:8],'b',.2)

    x=[1:8]';
    h=plot(x,mean(avg_acc_up,2),'r',x,mean(avg_hit_up,2),'g',x,mean(avg_fa_up,2),'b');
    legend('acc','hit','fa');
    for i=1:3
    set(h(i),'LineWidth',1.5);
    end
    ylim([-0.2 1]);
    title('averaged accuracy averaged across gender for upright faces','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    hold off
    
    subplot(3,1,2);
    ciplot(cilow_inv, cihigh_inv, [1:8],'r',.5)
    hold on
    ciplot(cilow_hit_inv, cihigh_hit_inv, [1:8],'g',.2)
    ciplot(cilow_fa_inv, cihigh_fa_inv, [1:8],'b',.2)

    x=[1:8]';
    h=plot(x,mean(avg_acc_inv,2),'r',x,mean(avg_hit_inv,2),'g',x,mean(avg_fa_inv,2),'b');
    legend('acc','hit','fa');
    for i=1:3
    set(h(i),'LineWidth',1.5);
    end
    ylim([-0.2 1]);
    title('averaged accuracy averaged across gender for inverted faces','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    hold off
    
    subplot(3,1,3);
    ciplot(cilow_diff, cihigh_diff, [1:8],'r',.5)
    hold on
    ciplot(cilow_hit_diff, cihigh_hit_diff, [1:8],'g',.2)
    ciplot(cilow_fa_diff, cihigh_fa_diff, [1:8],'b',.2)

    x=[1:8]';
    h=plot(x,mean(avg_acc_diff,2),'r',x,mean(avg_hit_diff,2),'g',x,mean(avg_fa_diff,2),'b');
    legend('acc','hit','fa');
    for i=1:3
    set(h(i),'LineWidth',1.5);
    end
    ylim([-0.2 1]);
    title('inversion effect across frequencies','Fontsize',18,'FontWeight', 'BOLD');
    set(gca,'XtickLabel',FreqLabels, 'Fontsize',12,'box','off')
    hold off


