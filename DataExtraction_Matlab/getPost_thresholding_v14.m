%% getSeizures
% gets post-ictal segments from iEEG portal for visualization
% NB: at the moment only works for segments of less than 2000 s
%TESTING!

function [] =  getPost_thresholding_v14(iPt)
%% Options and parameters
tic
% Parameters
sz_to_use = 10;             % number os seizures to use in this analysis (takes the first found)
search_times = [30 30 10 10 10 60 10 10 10 10 20 10 10 10 40];
search_time = search_times(iPt);       % min, length of post-ictal to analyse
low_f = 10;                 % Hz, lowest freq to consider in thresholding
high_f = 30;                % Hz, Highest Freq to consider in threhoslding
inter_length = 10;          % min, length of interictal period (for averaging purposes)
%              1    2    3   4     5    6    7     8   9    10   11
threshold = [4000 3500 4000 2000 10000 800 10000 2000 8000 4000 3000 10000 3000 2000 3000];   % thresholds based on eyballing of EEG traces 
sz_len_threshold = 25;      % s, cut-off between long and short seizures, determined by eye
F_Ord = 2;                  % Filter order
smooth_window_size = 5;     % s, length of the smoothing filter applied
median_window_size = 5;     % s, length of median filter

Fs = 400;                   % rounded Fs, for filtering
iCh = 1:16;                 
Type3 = 0;                  % Get type 3 seizures = 1

% ISI limitations
minISI = 5*60*60;           % s, Minimum ISI length, set to 5 Hr for tail seizures
maxISI = Inf*60*60;           % s, Maximum ISI length, set to inf for tail seizures

% training data time cutoff (days)
start_cutoff = 15*7;
end_cutoff = inf*7;

% Options
tails = true;              % enable to use tail seizures, otherwise lead will be used. More genrally, enable to select seizures based on time until next seizure
all_sz = true;             % use all valid seizures, if false, sz_to_use determines limit
plot_ind_sz = false;        % plot trace of seizures individually, best to use debugging so that one can show at a time
plot_scatter = false;       % plot sz_length vs post_length scatterplot
split_scatter = false;      % plot scatters for short and long seizures separately as well
plot_histogram = false;     % plot histogram of post-ictal lengths
plot_boxplot = false;       % plot boxplot for post-ictal times of long vs short seizures
plot_spectrogram = false;   % plot spectrograms for each channel of each seizure, recommend using debugging to view one seizure at a time
plot_all_on_one = false;    % plot all eeg traces on the one graph 
collect_inter = false;      % enable look at 1 hour of interictal instead, displaying mean and variance
threshold_average = true;  % Use the average interictal energy as the threshold
exp_window = false;         % Use an exponential smoothing window  
median_filter = true;      % applies a second smoothing window using a median filter
survival_curve = false;     % output a heatmap of post-ictal power ordered by seizure length
order_by_ISI = false;       % makes the 'survival' curve ordered by following ISI insead of seizure length, scatter becomes post_length vs ISI
ISI_compare = false;        % makes a scatter of sz vs post-ictal color coded for ISI length
next_sz = false;            % compares post_ictal to the length of the next seizure, only works if ISI is set to 0 to inf       

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

%% load information
Fs_actual = patient.data.sampleRate;
% fprintf('Actual freq = %d', Fs_actual)

load(['Portal Annots/' curPt '_Annots']);
load('Portal Annots/portalT0');

% chron. order
[SzTimes,I] = sort(SzTimes);
SzType = SzType(I);
SzDur = SzDur(I);


%% Select seizures to use
ISI = diff(SzTimes)/1e6;    %s, length of interseizure interval, ie time until next seizure

if tails
    ISI = [ISI minISI+1];  % Use tail seizures only
else
    ISI = [minISI+1 ISI]; % Use lead seizures only
end

%Remove type 3 seizures if not usung them
if ~Type3
    remove = SzType == 3;
    ISI(remove) = [];
    SzTimes(remove) = [];
    SzDur(remove) = [];
end

%Finds only lead/tail seizures in training period
SzDay = ceil(SzTimes/1e6/60/60/24);
training = SzDay > start_cutoff & SzDay < end_cutoff;
SzInd = find(ISI > minISI & ISI < maxISI & training);
ISI = ISI(SzInd);

%% Cycle through valid seizures
N = length(SzInd);

fprintf('\n%d seizures\n',N)

post_lengths = zeros(N,1);
sz_lengths = zeros(N,1);
averages = zeros(N,1);
stds = zeros(N,1);
power_traces = [];

if all_sz
    sz_to_use = N;
end


%%
%% LOOP STARTS HERE
%%

for n = 1:sz_to_use

    sz_lengths(n) = SzDur(SzInd(n));
    segment_length = search_time*60 + SzDur(SzInd(n));      % sec, Length of time to record
%     fprintf('\nSegment length = %d\n', segment_length)
    fprintf('Seizure %d of %d\n',n,N)
    % intialize start time
    t0 = SzTimes(SzInd(n))/1e6;                             % s, time to start recording

    
%% Collect data
    
    if collect_inter
        t0 = t0 + search_time*60;             % s, time to start recording interictal
        segment_length = inter_length*60;    % s, length of interictal from which to determine average
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
            power_traces = [power_traces; zeros(1, round(search_time*Fs_actual*60))];
            post_lengths(n) = inf;
            continue;
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
    
%% Extract 10-30Hz energy
    freq_range = [low_f high_f];                % HZ, Frequency range

    Wn = freq_range/(Fs/2);                     % Normalized cutoff
    [b, a] = butter(F_Ord, Wn, 'bandpass');     % define LP filter
    B_data = filtfilt(b,a,ZM_data);             % banded data looking at only key frequencies
    SQ_B_data = B_data.^2;                      % Convert to energy
    
%% Filters/Smoothing 
    
    if exp_window
        t_half = 1/4;                                           % fraction of smoothing window, half life of exponential window's decay
        alpha = 1/t_half * log(2) / (smooth_window_size);       % exponential coefficient
        t_w = 1/Fs:1/Fs:smooth_window_size;                     % s, timescale for the window
        b = exp(-alpha * t_w);                                  % exponential window
        a = sum(b);
        F_SQ_B_data = filter(b,a,SQ_B_data);                    % Moving average of squared data => power.
        
        P_data = F_SQ_B_data(Fs * smooth_window_size:end, :);   % shift the output of the filter 20s (ie window size) to the left, now the ouput of the filter is on its leftmost input
     
    else
        b = ones(2*Fs*smooth_window_size,1);
%         b = hamming(2*Fs * smooth_window_size);             % filter variables
        b = b(Fs * smooth_window_size+1:end);               % just take the right hand side of the window
        a = sum(b);
        P_data = filter(b,a,SQ_B_data);                     % Moving average of squared data => power.
        P_data = P_data(Fs * smooth_window_size:end, :);    % adjusts for using a one-sided filter
    end
    
    T_P_data = sum(P_data,2);                       % Sum across channels to get total power
    
    if median_filter
        T_P_data = medfilt1(T_P_data,Fs*median_window_size);
    end
    
    high_ix = find(T_P_data>10000);                 % Find indicies where power is blown out
    if ~collect_inter
        T_P_data(high_ix) = 10000;                  % Cap power at 10000, useful for visualization
        if size(T_P_data,1)< round(search_time*Fs_actual*60)
            T_P_data(round(search_time*Fs_actual*60)) = 0;      % adds zeros at the end of short T_P_data
        end
        power_traces = [power_traces; T_P_data(max(1,end-(round(search_time*Fs_actual*60)-1)):end)']; %this causes matrix dimention mismatch when T_P_data is too short
    end
    
    sz_end = SzDur(SzInd(n));                                   % s, endpoint of seizure
    post_start = round(sz_end*Fs);    % timebins, point to start looking for end of post
    
    
    
%% Generate post-ictal lengths based on arbitrary thresholding of power


    if threshold_average
        av_str = ['thresholds/averages_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600)];
        load(av_str);                   % loaded matrix should be called 'averages'
        try
            ix = find(T_P_data(post_start:end)>averages(n), 1);        % Find first index past threhsold
        catch
            ix = nan;
        end

    else
%         ix = find(T_P_data(post_start:end)>threshold(iPt), 1);        % Find first index past threhsold
        
%         idx=T_P_data(post_start:end)>threshold(iPt);
    end
    % UNUSED code for continuous values above threshold
%     ii1=strfind([0 idx' 0],[0 1]);
%     ii2=strfind([0 idx' 0],[1 0]);
%    	ii=(ii2-ii1)>=Fs*10;
%    % out=arrayfun(@(x,y) A(x:y),ii1(ii),ii2(ii),'un',0);
%    % celldisp(out)
%    
%   	try
%      	ix = ii1(find(ii),1);
%     catch
% %     	continue;
%     end 
    
    try
        post_lengths(n) = ix/Fs;       % s, time till above threshold
    catch
        post_lengths(n) = inf;         % use inf as a flag to delete this sz from list
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
            t_plot = t(1:end - smooth_window_size*Fs+1);
            plot(t_plot, T_P_data, 'r')                      % Overlay power values, shown in red
            
            plot([sz_end sz_end], [0 8000])             % Plot seizure end time
            
            post_end = post_lengths(n) + sz_end;        % s, time postictal ends relative to start of seizure
            
            plot([post_end post_end],[0 8000], 'g')     % Plot end of post-ictal
            plot([0 search_time*60],[(averages(n)) (averages(n))], 'g')
            hold off
            title(['Pt ' num2str(iPt) ', exp = ' num2str(smooth_window_size) ', med_len = ' num2str(median_window_size)])
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
%% Remove inf length post-ictals

ix = find(post_lengths~=inf);             %Sz indecies for which post-ictal length > 0
post_lengths = post_lengths(ix);        % Remove zero-length post-ctal times = ['Post_Ictal_Lengths/Pt_' num2str(iPt)]
savepath = ['Post_Ictal_Lengths/Pt_' num2str(iPt)];
save(savepath, 'post_lengths');             % Sav post_ictal lengths

R_sz_lengths=sz_lengths(ix);              % Remove seizures of inf length post-ictal times
R_ISI = ISI(ix);                            % Remove ISI lengths for invalid post-ictal
sz_shown = size(post_lengths,1);            % Number of seizures used (excluding ones with zero length post-ictal

fprintf([num2str(sz_shown) ' of ' num2str(N) ' seizures have valid post-ictal times\n'])

%% Histogram
if plot_histogram
    figure
    histogram(post_lengths,20);             % Distribution of post-ictal lengths
    title(['Pt - ' num2str(iPt)])
    xlabel('Post-ictal length (s)')
end

%% Scatter plot for all seizures
if plot_scatter
    figure(1);
    if order_by_ISI
        scatter(post_lengths, R_ISI/3600);      % Scatterplot of seizure length vs post-ictal length.
        axis([0 max(post_lengths) 0 max(R_ISI/3600)]);
        lsline;                                 % adds line of best fit
        title(['Pt - ' num2str(iPt)], 'Fontname', 'Calibri', 'Fontsize', 10);
        xlabel('Post-ictal Length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
        ylabel('ISI length (hr)', 'Fontname', 'Calibri', 'Fontsize', 6);
        set(gca, 'FontName', 'Calibri', 'Fontsize', 6);

        [R,P] = corrcoef(R_ISI/3600, post_lengths);   % Correlation coefficient and associated p-value
        %Display R-sqr and P to console
        fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
        fprintf('p=value = %d\n', P(1,2));
    %     savefig(1,['Figures/ScatSurv/scatter_' num2str(iPt)]);
        fig=gcf;    % get current figure)
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 3 3];    % Image size in inches
        print(['Figures/ScatSurv/scatter_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600) '_ISI'], '-dtiff', '-r600')  % -rx changes dpi/ppi to x
    else
        if next_sz
            R_sz_lengths = R_sz_lengths(2:end);
            post_lengths = post_lengths(1:end-1);
        end
        scatter(R_sz_lengths, post_lengths);      % Scatterplot of seizure length vs post-ictal length.
        axis([0 max(R_sz_lengths) 0 max(post_lengths)])
        lsline;                                 % adds line of best fit
        title(['Pt - ' num2str(iPt)], 'Fontname', 'Calibri', 'Fontsize', 10);
        xlabel('Seizure Length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
        ylabel('Post-ictal length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
        set(gca, 'FontName', 'Calibri', 'Fontsize', 6);

        [R,P] = corrcoef(post_lengths,R_sz_lengths);   % Correlation coefficient and associated p-value
        %Display R-sqr and P to console
        fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
        fprintf('p=value = %d\n', P(1,2));
    %     savefig(1,['Figures/ScatSurv/scatter_' num2str(iPt)]);
        fig=gcf;    % get current figure)
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 3 3];    % Image size in inches
        if next_sz
            print(['Figures/ScatSurv/scatter_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600) '_next_sz'], '-dtiff', '-r600')  % -rx changes dpi/ppi to x
        else
            print(['Figures/ScatSurv/scatter_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600)], '-dtiff', '-r600')  % -rx changes dpi/ppi to x
        end
    end
    
    if split_scatter
    %%  Short seizures scatter
        figure;
        ix = find(R_sz_lengths<sz_len_threshold);
        scatter(R_sz_lengths(ix), post_lengths(ix));      % Scatterplot of seizure length vs post-ictal length.
        axis([0 sz_len_threshold 0 max(post_lengths(ix))])
        lsline;                                 % adds line of best fit
        title(['Pt - ' num2str(iPt)])
        xlabel('Seizure Length (s)')
        ylabel('Post-ictal length (s)')

        [R,P] = corrcoef(post_lengths(ix),R_sz_lengths(ix));   % Correlation coefficient and associated p-value
        %Display R-sqr and P to console
        fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
        fprintf('p=value = %d\n', P(1,2));

    %%  Long Seizures only
        figure;
        ix = find(R_sz_lengths>=sz_len_threshold);
        scatter(R_sz_lengths(ix), post_lengths(ix));      % Scatterplot of seizure length vs post-ictal length.
        axis([sz_len_threshold max(R_sz_lengths(ix)) 0 max(post_lengths)])
        lsline;                                 % adds line of best fit
        title(['Pt - ' num2str(iPt)])
        xlabel('Seizure Length (s)')
        ylabel('Post-ictal length (s)')

        [R,P] = corrcoef(post_lengths(ix),R_sz_lengths(ix));   % Correlation coefficient and associated p-value
        %Display R-sqr and P to console
        fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
        fprintf('p=value = %d\n', P(1,2));
    end
end

%% Create scatter with 

if ISI_compare

    short = find(R_ISI<2*3600);
    medium = find(R_ISI<5*3600 & R_ISI>2*3600);
    long = find(R_ISI>5*3600);
    sz_short = R_sz_lengths(short);
    sz_med = R_sz_lengths(medium);
    sz_long = R_sz_lengths(long);
    
    post_short = post_lengths(short);
    post_med = post_lengths(medium);
    post_long = post_lengths(long);
    
    figure
    scatter(sz_long, post_long, 'r');
    hold on;
    scatter(sz_med, post_med, 'k');
    scatter(sz_short, post_short, 'b');
    hold off;
    
    axis([0 max(R_sz_lengths) 0 max(post_lengths)])
    title(['Pt - ' num2str(iPt)], 'Fontname', 'Calibri', 'Fontsize', 10);
    xlabel('Seizure Length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
    ylabel('Post-ictal length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
    set(gca, 'FontName', 'Calibri', 'Fontsize', 6);
    legend('Long', 'Intermediate', 'Short')

    [R,P] = corrcoef(post_lengths,R_sz_lengths);   % Correlation coefficient and associated p-value
    %Display R-sqr and P to console
    fprintf('\nR-square = %d\n', R(1,2)*R(1,2));
    fprintf('p=value = %d\n', P(1,2));
    %     savefig(1,['Figures/ScatSurv/scatter_' num2str(iPt)]);
    fig=gcf;    % get current figure)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 3 3];    % Image size in inches
    print(['Figures/ScatSurv/scatter_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600) '_3color'], '-dtiff', '-r600')  % -rx changes dpi/ppi to x


end

%% Box Plot
if plot_boxplot
    sz_time_threshold = 25;                 % s, cutoff time between long and short seizures

    short_sz = post_lengths(R_sz_lengths<sz_time_threshold);   % post ictal lengths of short seizures
    long_sz = post_lengths(R_sz_lengths>=sz_time_threshold);   % post ictal lengths of long seizures

    figure
    grp = [zeros(size(short_sz')), ones(size(long_sz'))];           % defines the size of the datasets
    boxplot([short_sz',long_sz'], grp, 'Labels',{'Short','Long'})   % plot the boxplots
    title(['Patient ' num2str(iPt)]);
    
    % t test
    [h,p,ci] = ttest2(short_sz, long_sz) % Test for significant difference bw post-ictal populations
end

%% Survival Curve
if survival_curve
    if threshold_average
        S_power_traces = power_traces./averages(1:sz_to_use);
    else
        S_power_traces = power_traces/threshold(iPt);
    end
    

    % grab and order valid post-ictals
    S_power_traces(isnan(S_power_traces))=0;
    R_S_p_t = S_power_traces(find(sum(S_power_traces,2)~=0),:);     % removes sz with excess dropout (previous set to all 0)
    
    
    
    if order_by_ISI
        R_ISI = ISI(find(sum(S_power_traces,2)~=0));                % removes sz with excess dropout
        [O_R_ISI, ix] = sort(R_ISI');
    else
        R_sz_lengths = sz_lengths(find(sum(S_power_traces,2)~=0),:);    % removes sz with excess dropout (previous set to all 0)
        [O_R_sz_lengths,ix] = sort(R_sz_lengths);                       % get indicies of seizure sorted by length
    end
    O_R_S_p_t = R_S_p_t(ix,:);                                      % sort post-ictal using seizure (or ISI) length  (ascending)
    
    O_R_S_p_t(find(O_R_S_p_t>2)) = 2;                           % limit intensity to twice average
    
    
    figure(2);
    
    %Plot sz lengths
    subplot(1,2,1)                          % plot to left of figure
    if order_by_ISI
        sz_num = 1:size(O_R_ISI,1);                              % make y-axis values (sz number) 
        scatter(O_R_ISI, sz_num, 3, 'filled')             % plot sz_lengths
        axis([0 max(O_R_ISI)*1.1 .5 size(O_R_ISI,1)+.5]); % adjust axis to fit best
    else
        sz_num = 1:size(O_R_sz_lengths,1);                              % make y-axis values (sz number) 
        scatter(O_R_sz_lengths, sz_num, 3, 'filled')             % plot sz_lengths
        axis([0 max(O_R_sz_lengths)*1.1 .5 size(O_R_sz_lengths,1)+.5]); % adjust axis to fit best
    end
    set(gca, 'ycolor', 'none', 'ydir', 'reverse', 'position', [0.03 0.1 0.15 0.8]); % ycolor removes y axes, yir and xdir reverses axes, position sets position in [left, bottom, width, height]
    xlabel('Seizure length (s)', 'Fontname', 'Calibri', 'Fontsize', 6);
    set(gca, 'FontName', 'Calibri', 'Fontsize', 4);
    
    % Plot 'survival' curve with blue -> red colorbar
    subplot(1,2,2)      % allocate to right of image
    if order_by_ISI
        imagesc([1/(Fs_actual*60) size(O_R_S_p_t, 2)/(Fs_actual*60)], [1 size(O_R_ISI, 2)], O_R_S_p_t);
    else
        imagesc([1/(Fs_actual*60) size(O_R_S_p_t, 2)/(Fs_actual*60)], [1 size(O_R_sz_lengths, 2)], O_R_S_p_t);
    end
    set(gca, 'ycolor', 'none', 'box', 'off', 'position', [0.20 0.1 0.68 0.8]); % y color to none removes axis, box off removes top x axis, position sets position in [left, bottom, width, height]
    xlabel('Time since seizure end (min)', 'Fontname', 'Calibri', 'Fontsize', 6);
    set(gca, 'FontName', 'Calibri', 'Fontsize', 4) % gca = get current axis
    
    % Define color map
    cmap = zeros(41,3);
    cmap(:,1) = [0:      1/20:         1   1-92/(20*256):  -92/(20*256):  164/256]';          % Red   fades up to white, then solid
    cmap(:,2) = [34/256: 222/(20*256): 1   1-253/(20*256): -253/(20*256): 3/256]';          % Green  fades jup to white then fadces down
    cmap(:,3) = [88/256: 168/(20*256): 1   1-203/(20*256): -203/(20*256): 53/256]';      	% Blue   solid then fades down
    colormap(cmap);                             % sets colormap to above
    colorbar('Position', [.9 0.096 .04 .76]);    % shows colormap and sets it s position
    if order_by_ISI
        h = suptitle(['Pt ' num2str(iPt) ' ISI length, ISI from ' num2str(minISI/3600) ' to ' num2str(maxISI/3600) ' hours']);         % Set title for whole graph
    else
        h = suptitle(['Pt ' num2str(iPt) ' seizure length, ISI from ' num2str(minISI/3600) ' to ' num2str(maxISI/3600) ' hours']);         % Set title for whole graph
    end
    set(h,'Fontsize',8,'Fontname','Calibri')    % Alter title font
    % savefig(2, ['Figures/ScatSurv/surv_' num2str(iPt)]);
    % Save fig as tif with set resolution and size
    fig=gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 3 2];    % Image size in inches
    if order_by_ISI
        print(['Figures/ScatSurv/surv_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600) '_byISI'], '-dtiff', '-r600')  % -rx changes dpi/ppi to x
    else
        print(['Figures/ScatSurv/surv_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600)], '-dtiff', '-r600')  % -rx changes dpi/ppi to x
    end

end


if collect_inter
    filename = ['thresholds/averages_' num2str(iPt) '_' num2str(minISI/3600) '_' num2str(maxISI/3600)];
    save(filename, 'averages');
end

x=1;    % just for debugging, can delete
toc
