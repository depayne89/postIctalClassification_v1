%% getSeizures
% gets post-ictal segments from iEEG portal for visualization
% NB: at the moment only works for segments of less than 2000 s
%TESTING!

function [] =  getPost_thresholding_v6(iPt)

% IEEG LOGIN HERE
login = 'depayne';
pword = 'pwd.bin';

% Patients
Patient{1} = '23_002';
Patient{2} = '23_003';
Patient{3} = '23_004';
Patient{4} = '23_005';
Patient{5} = '23_006';
Patient{6} = '23_007';

Patient{7} = '24_001';
Patient{8} = '24_002';
Patient{9} = '24_004';
Patient{10} = '24_005';

Patient{11} = '25_001';
Patient{12} = '25_002';
Patient{13} = '25_003';
Patient{14} = '25_004';
Patient{15} = '25_005';

temp = ['Sz_' num2str(iPt)];
parent_path = 'C:/Users/depayne/Desktop/NVData/Pt13_09_01_17/';
mkdir(parent_path,temp);
save_path = [temp '/'];

%% link to portal session
curPt = Patient{iPt};
patient = IEEGSession(['NVC1001_' curPt '_2'],login,pword);

%% parameters
Fs_actual = patient.data.sampleRate;
fprintf('Actual freq = %d', Fs_actual)
Fs = 400;  % for filtering
iCh = 1:16;

% Only save seizures with at least XX as seizure-free time beforehand (lead
% time)
% Lead time to include the seizures in seconds
leadTime = 5*60*60;

% Get type 3 seizures = 1
Type3 = 0;

% training data time cutoff (days)
start_cutoff = 15*7;
end_cutoff = inf*7;

%% load information
load(['Portal Annots/' curPt '_Annots']);
load('Portal Annots/portalT0');

% chron. order
[SzTimes,I] = sort(SzTimes);
SzType = SzType(I);
SzDur = SzDur(I);


%% Select only tail seizures after 100 days
ISI = diff(SzTimes)/1e6;    %s, length of interseizure interval, ie time until next seizure
ISI = [leadTime+1 ISI];

%Remove type 3 seizures if not usung them
if ~Type3
    remove = SzType == 3;
    ISI(remove) = [];
    SzTimes(remove) = [];
    SzDur(remove) = [];
end

%Finds only lead seizures in training period
SzDay = ceil(SzTimes/1e6/60/60/24);
training = SzDay > start_cutoff & SzDay < end_cutoff;
SzInd = find(ISI > leadTime & training);

%% Cycle through valid seizures
N = length(SzInd);

fprintf('\n%d seizures\n',N)

post_lengths = zeros(N,1);
sz_lengths = zeros(N,1);

%For every seizure
for n = 1:N
    sz_lengths(n) = SzDur(SzInd(n));
    search_time = 5;                                    % min, length of post-ictal time to look at
    segment_length = search_time*60 + SzDur(SzInd(n));   % sec, Length of time to record
    fprintf('\nSegment length = %d\n', segment_length)
    fprintf('Seizure %d of %d\n',n,N)
    % intialize start time
    t0 = SzTimes(SzInd(n))/1e6;                          % sec, time to start recording
    
%% Collect data
    try
        Data = getvalues(patient.data,t0 * 1e6,segment_length * 1e6,iCh);
    catch
        % try again
        try
            Data = getvalues(patient.data,t0 * 1e6,segment_length * 1e6,iCh);
        catch
            % maybe lost connection.
            display('early termination');
            continue;
        end
    end
    
%% Remove seziures with excessive dropout
    %pre-filter
    Data(isnan(Data(:,1)),:) = 0;                   % set all nan to zero
    zero_frac = numel(find(Data==0))/numel(Data);    % fraction of data than = 0
    if zero_frac > 0.1
       fprintf('\n Excess Dropout, skipping this seizure');
       continue;
    end

%% Filter, zero-mean and offset traces
    Fs=400;                     % Approx sample rate of EEG
    NChs = 16;                  % Number of channels
    NSamples = size(Data,1);     % Number of samples
    
    t = 1/Fs:1/Fs:NSamples/Fs;  % Time vector in seconds
    
    %Zero mean the data instead of highpass filter
    ZM_data = Data - repmat(mean(Data,1)',1,NSamples)';
    
    F_Ord = 2;                  % Filter Order
    Fc = 120;                   % Hz, filter cut off
    Wn = Fc/(Fs/2);             % Normalized cutoff
    [b, a] = butter(F_Ord, Wn); % define LP filter
    
    F_ZM_data = filtfilt(b,a,ZM_data);  % Apply lowpass filter
    
    Ch_offset = 400;                                    % Distance between channels
    Offset_vector = Ch_offset:Ch_offset:Ch_offset*NChs; % Offset for each channel
    Offset_Mat = repmat(Offset_vector', 1, NSamples);   % Extend vector to a matrix
    
    OS_F_ZM_data = F_ZM_data + Offset_Mat';             % Add offset to the data
    
%% Extract 10-30Hz power
    low_f = 10;                                 % Hz, lowest freq to consider in thresholding
    high_f = 20;                                % Hz, Highest Freq to consider in threhoslding
    freq_range = [low_f high_f];                % HZ, Frequency range
    
    F_Ord = 2;                                  % Filter Order
    Wn = freq_range/(Fs/2);                     % Normalized cutoff
    [b, a] = butter(F_Ord, Wn, 'bandpass');     % define LP filter
    B_data = filtfilt(b,a,ZM_data);             % banded data looking at only key frequencies
    
%% Calculate power
    SQ_B_data = B_data.^2;                      % Convert to energy

    window_size = round(Fs*20);                   % samples, time power is calculated over
    b = (1/window_size) * ones(1,window_size);  % filter variables
    a = 1;

    P_data = filter(b,a,SQ_B_data);          % Moving average of squared data => power. 
    T_P_data = sum(P_data,2);                % Sum across channels to get total power
    high_ix = find(T_P_data>10000);          % Find indicies where power is blown out
    T_P_data(high_ix) = 10000;               % Cap power at 10000
    
%% Generate post-ictal lengths based on arbitrary thresholding of power


    sz_end = SzDur(SzInd(n));                                   % s, endpoint of seizure
    post_start = sz_end*Fs + 1.2*window_size;                   % timebins, point to start looking for end of post
    %              1    2    3   4     5    6    7     8   9    10   11
    threshold = [4000 3500 4000 2000 10000 800 10000 2000 8000 4000 3000 10000 3000 2000 3000];   % thresholds based on eyballing of EEG traces 
    ix = find(T_P_data(post_start:end)>threshold(iPt), 1);        % Find first index past threhsold

    try
        post_lengths(n) = (ix + 1.2*window_size)/Fs;       % s, time till above threshold
    catch
        post_lengths(n) = 0;
    end    
    
%% Plot EEG traces with power overlaid 
%     figure;
%     sz_end = SzDur(SzInd(n))                    % s, endpoint of seizure
%     plot(t,OS_F_ZM_data, 'k');                  % Plot EEG channels
%     hold on
%     plot(t, T_P_data, 'r')                      % Overlay power values, shown in red
%     plot([sz_end sz_end], [0 8000])             % Plot seizure end time
%     post_end = post_lengths(n) + sz_end;        % s, time postictal ends relative to start of seizure
%     plot([post_end post_end],[0 8000], 'g')     % Plot end of post-ictal
%     hold off

 end % end seizure loop
%% Remove zero length post-ictals

ix = find(post_lengths~=0);             %Sz indecies for which post-ictal length > 0
post_lengths = post_lengths(ix);        % Remove zero-length post-ctal times
sz_lengths=sz_lengths(ix);              % Remove seizures of zero length post-ictal times


%% Plot Figures

[R,P] = corrcoef(post_lengths,sz_lengths);   % Correlation coefficient and associated p-value

%Display R-sqr and P to console
fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
fprintf('p=value = %d\n', P(1,2));

sz_shown = size(post_lengths,1);            % Number of seizures used (excluding ones with zero length post-ictal

fprintf([num2str(sz_shown) ' of ' num2str(N) ' seizures have valid post-ictal times\n'])

figure
histogram(post_lengths,20);             % Distribution of post-ictal lengths
title(['Pt - ' num2str(iPt)])
xlabel('Post-ictal length (s)')

figure
scatter(sz_lengths, post_lengths);      % Scatterplot of seizure length vs post-ictal length.
lsline;                                 % adds line of best fit
title(['Pt - ' num2str(iPt)])
xlabel('Seizure Length (s)')
ylabel('Post-ictal length (s)')
