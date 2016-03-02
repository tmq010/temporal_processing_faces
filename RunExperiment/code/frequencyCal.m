clear all
close all

%% to calculate the actual presentation frequency for 8 hertz stimuli
time8=[1/8:1/8:1]; % ideal 8 hertz
timeRef=[1/60:1/60:1]; % all possible presnetation time for a monitor that has refreshing rate as 60 hz

for i = 1:8
tmp = abs(time8(i)-timeRef);
[x x] = min(tmp);
timeAct(i) = timeRef(x); % timeAct is the actual presentation time (s) that you get 
timeIdx(i) = x;
end

%% what are the possible frequencies for a 60hz monitor?

for i=1:60
    possibleFreq(i,1) = 60/i;
    possibleFreq(i,2) = i;
end


%% sinusoid binned by different numbers of frames

clx
t =[0:1:40]; % Time Samples
f =500; % Input Signal Frequency
fs = 500*20; % Sampling Frequency
x = sin(2*pi*f/fs*t); % Generate Sine Wave  
figure(1);
stem(t,x,'r'); % View the samples
figure(2);
stem(t*1/fs*1000,x,'r'); % View the samples
hold on;
plot(t*1/fs*1000,x); % Plot Sine Wave



