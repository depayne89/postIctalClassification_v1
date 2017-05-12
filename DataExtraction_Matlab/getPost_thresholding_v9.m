%% getSeizures
% gets post-ictal segments from iEEG portal for visualization
% NB: at the moment only works for segments of less than 2000 s
%TESTING!

function [] =  getPost_thresholding_v9(iPt)
%% Options and parameters

% Parameters
sz_to_use = 10;          % number os seizures to use in this analysis (takes the first found)
search_time = 5;       % min, length of post-ictal to analyse
low_f = 10;             % Hz, lowest freq to consider in thresholding
high_f = 30;            % Hz, Highest Freq to consider in threhoslding
inter_start = 60;        % min, time until start of interictal period (for averaging purposes)
%              1    2    3   4     5    6    7     8   9    10   11
threshold = [4000 3500 4000 2000 10000 800 10000 2000 8000 4000 3000 10000 3000 2000 3000];   % thresholds based on eyballing of EEG traces 
sz_len_threshold = 25;  % s, cut-off between long and short seizures, determined by eye
F_Ord = 2;              % Filter order

% Options
tails = true;              % enable to use tail seizures, otherwise lead will be used
all_sz = true;             % use all valid seizures, if false, sz_to_use determines limit
plot_ind_sz = false;        % plot trace of seizures individually, best to use debugging so that one can show at a time
plot_scatter = true;       % plot sz_length vs post_length scatterplot
plot_histogram = false;     % plot histogram of post-ictal lengths
plot_boxplot = false;       % plot boxplot for post-ictal times of long vs short seizures
plot_spectrogram = false;   % plot spectrograms for each channel of each seizure, recommend using debugging to view one seizure at a time
plot_all_on_one = false;    % plot all eeg traces on the one graph 
collect_inter = false;      % enable look at 1 hour of interictal instead, displaying mean and variance
threshold_average = true;  % Use the average interictal energy as the threshold

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


%% Select seizures to use
ISI = diff(SzTimes)/1e6;    %s, length of interseizure interval, ie time until next seizure

if tails
    ISI = [ISI leadTime+1];  % Use tail seizures only
else
    ISI = [leadTime+1 ISI]; % Use lead seizures only
end

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
averages = zeros(N,1);
stds = zeros(N,1);


if all_sz
    sz_to_use = N;
end

%For every seizure
for n = 1:sz_to_use
    sz_lengths(n) = SzDur(SzInd(n));
    segment_length = search_time*60 + SzDur(SzInd(n));      % sec, Length of time to record
%     fprintf('\nSegment length = %d\n', segment_length)
    fprintf('Seizure %d of %d\n',n,N)
    % intialize start time
    t0 = SzTimes(SzInd(n))/1e6;                             % s, time to start recording

    
%% Collect data
    
    if collect_inter
        t0 = t0 + inter_start*60;             % s, time to start recording interictal
        segment_length = search_time*60;    % s, length of interictal from which to determine average
    end
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
       if collect_inter
            Data = Data(sum(Data,2)>1e-10,:);
            fprintf('\nExcess Dropout\n');
       else
            fprintf('\n Excess Dropout, skipping this seizure\n');
%             continue;
       end
    end
    if size(Data,1)==0
        continue
    end

%% Filter, zero-mean and offset traces
    
    Fs=400;                     % Approx sample rate of EEG
    NChs = 16;                  % Number of channels
    NSamples = size(Data,1);    % Number of samples
    
    t = 1/Fs:1/Fs:NSamples/Fs;  % Time vector in seconds
    
    %Zero mean the data instead of highpass filter
    ZM_data = Data - repmat(mean(Data,1)',1,NSamples)';
    
    Fc = 120;                   % Hz, filter cut off
    Wn = Fc/(Fs/2);             % Normalized cutoff
    [b, a] = butter(F_Ord, Wn); % define LP filter
    
    F_ZM_data = filtfilt(b,a,ZM_data);  % Apply lowpass filter
    
    Ch_offset = 400;                                    % Distance between channels
    Offset_vector = Ch_offset:Ch_offset:Ch_offset*NChs; % Offset for each channel
    Offset_Mat = repmat(Offset_vector', 1, NSamples);   % Extend vector to a matrix
    
    OS_F_ZM_data = F_ZM_data + Offset_Mat';             % Add offset to the data
    
%% Extract 10-30Hz power
    freq_range = [low_f high_f];                % HZ, Frequency range

    Wn = freq_range/(Fs/2);                     % Normalized cutoff
    [b, a] = butter(F_Ord, Wn, 'bandpass');     % define LP filter
    B_data = filtfilt(b,a,ZM_data);             % banded data looking at only key frequencies
    
%% Calculate power
    SQ_B_data = B_data.^2;                   % Convert to energy

    window_size = round(Fs*20);                    % samples, time power is calculated over
    
    b = hamming(window_size)*1;              % filter variables
    a = sum(hamming(window_size));
    P_data = filter(b,a,SQ_B_data);           % Moving average of squared data => power. 
    T_P_data = sum(P_data,2);                   % Sum across channels to get total power
    high_ix = find(T_P_data>10000);             % Find indicies where power is blown out
    if ~collect_inter
        T_P_data(high_ix) = 10000;                  % Cap power at 10000, useful for visualization
    end
%% Generate post-ictal lengths based on arbitrary thresholding of power


    sz_end = SzDur(SzInd(n));                                   % s, endpoint of seizure
    post_start = round(sz_end*Fs + window_size);                   % timebins, point to start looking for end of post
    
    if threshold_average
        av_str = ['averages_' num2str(iPt)];
        load(av_str);
        av_av = mean(averages);
        ix = find(T_P_data(post_start:end)>av_av, 1);
    else
        ix = find(T_P_data(post_start:end)>threshold(iPt), 1);        % Find first index past threhsold
    end
    try
        post_lengths(n) = (ix + 1.2*window_size)/Fs;       % s, time till above threshold
    catch
        post_lengths(n) = 0;
    end    
    
%% Plot individual seizures
    
    if plot_all_on_one
        sz_end = SzDur(SzInd(n));                    % s, endpoint of seizure   
        t_post = t(sz_end*Fs:end) - sz_end;         % adjusted timeseries from end of seizure
        T_P_post = T_P_data(sz_end*Fs:end);         % adjusted sample from end of seizure
        if sz_end< sz_len_threshold                 % display plots in different colors based on time threshold
            plot(t_post, T_P_post, 'k')             % Overlay power values, shown in red
        else
            plot(t_post, T_P_post, 'r')
        end
%         post_end = post_lengths(n);        % s, time postictal ends relative to start of seizure
        plot([0 search_time*60],[threshold(iPt) threshold(iPt)], 'g')     % Plot end of post-ictal
        hold on
        
    else
        
        if plot_ind_sz
            figure;
            sz_end = SzDur(SzInd(n));                    % s, endpoint of seizure
            plot(t,OS_F_ZM_data, 'k');                  % Plot EEG channels
            hold on
            plot(t, T_P_data, 'r')                      % Overlay power values, shown in red
            plot([sz_end sz_end], [0 8000])             % Plot seizure end time
            post_end = post_lengths(n) + sz_end;        % s, time postictal ends relative to start of seizure
            plot([post_end post_end],[0 8000], 'g')     % Plot end of post-ictal
            hold off
        end
    
        if plot_spectrogram
             figure;
            for ch = iCh

            %     t = .0025:.0025:30;
            %     length = size(t,2);
            %     y = rand(size(t)) +10 * sin(200 * pi * t).*rand(size(t));

                y = Data(:,ch);

                specSize = 32;                                      % size of output images
                nfft = ceil((4/3)*(-1+(1+1.5*NSamples)^0.5));       % derived to make a square spectrogram (see IBM notes)
                overlap = round(nfft/4);                            % overlap between timebins
                [S, F, T] = spectrogram(y,nfft,overlap,nfft, Fs);   %
                maxHz = 50;
                S = imresize(abs(S(1:ceil(size(S,1)*2*maxHz/Fs),:)), [specSize, specSize]);

                xvec = 0:search_time/32:search_time;
                yvec = 0:maxHz/32:maxHz;
                subplot(4,4,ch), imagesc(xvec,yvec,S);
        %         set(gca,'XtickLabel',[],'YtickLabel',[]);
            end
        end
    end
    
    if collect_inter
        ZR_Data = T_P_data(find(T_P_data>1e-10));
        averages(n) = mean(T_P_data);
        stds(n) = std(T_P_data);
%         length_remainig
    end
    
 end % end seizure loop
%% Remove zero length post-ictals

ix = find(post_lengths~=0);             %Sz indecies for which post-ictal length > 0
post_lengths = post_lengths(ix);        % Remove zero-length post-ctal times
sz_lengths=sz_lengths(ix);              % Remove seizures of zero length post-ictal times


%% Plot Figures

[R,P] = corrcoef(post_lengths,sz_lengths);   % Correlation coefficient and associated p-value

%Display R-sqr and P to console
if ~collect_inter
    fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
    fprintf('p=value = %d\n', P(1,2));
end

sz_shown = size(post_lengths,1);            % Number of seizures used (excluding ones with zero length post-ictal

fprintf([num2str(sz_shown) ' of ' num2str(N) ' seizures have valid post-ictal times\n'])

% Histogram
if plot_histogram
    figure
    histogram(post_lengths,20);             % Distribution of post-ictal lengths
    title(['Pt - ' num2str(iPt)])
    xlabel('Post-ictal length (s)')
end

% Scatter Plot
if plot_scatter
    figure
    scatter(sz_lengths, post_lengths);      % Scatterplot of seizure length vs post-ictal length.
    lsline;                                 % adds line of best fit
    title(['Pt - ' num2str(iPt)])
    xlabel('Seizure Length (s)')
    ylabel('Post-ictal length (s)')
end

% Box Plot
if plot_boxplot
    sz_time_threshold = 25;                 % s, cutoff time between long and short seizures

    short_sz = post_lengths(sz_lengths<sz_time_threshold);   % post ictal lengths of short seizures
    long_sz = post_lengths(sz_lengths>=sz_time_threshold);   % post ictal lengths of long seizures

    figure
    grp = [zeros(size(short_sz')), ones(size(long_sz'))];           % defines the size of the datasets
    boxplot([short_sz',long_sz'], grp, 'Labels',{'Short','Long'})   % plot the boxplots
    title(['Patient ' num2str(iPt)]);
    
    % t test
    [h,p,ci] = ttest2(short_sz, long_sz) % Test for significant difference bw post-ictal populations
end
x=1;    % just for debugging, can delete

