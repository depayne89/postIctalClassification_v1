

% post-ictal viewer

clc
clear
close all

% create an array of random numbers as a surrogate for real eeg
NChs = 16;                      % number of EEG channels
Fs = 400;                       % sample rate of EEG
T = 1000;                       % number of seconds in entire data range
NSamples = T*Fs;                % number of samples in entire data range
data = randn(NSamples,NChs);    % mock data (zero mean)

ZM_data = data - repmat(mean(data,1)',1,NSamples)';

F_Ord = 2;                      % filter order
Fc = 35;                        % Hz, filter cut off freq
Wn = Fc/(Fs/2);                 % normalised cutoff
[b, a] = butter(F_Ord, Wn);     % def LP filter

ZM_F_data = filtfilt(b,a,ZM_data);

Ch_Offset = 4;                  % offset for each channel
Offset_Vector = 0:Ch_Offset:Ch_Offset*(NChs-1);
Offset_Mat = repmat(Offset_Vector',1,NSamples)';

OS_ZM_F_data = ZM_F_data + Offset_Mat;

Plot_T = 20;                    % s, time for plot
Plot_Samples = Plot_T*Fs;
Plot_Advance = 5;               %s, time to advance figure
Samples_advance = Plot_Advance*Fs;

fig1 = figure;
ax1 = axes(fig1);
for n=1:Samples_advance:NSamples-Plot_Samples
    plot(ax1,OS_ZM_F_data(n:n+Plot_Samples,:),'k')
    axis off
    pause
end

