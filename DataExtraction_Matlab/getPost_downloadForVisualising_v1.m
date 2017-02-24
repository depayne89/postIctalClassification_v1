%% getSeizures
% gets pre-ictal segments from iEEG portal
% NB: at the moment only works for segments of less than 2000 s
%TESTING!

function [] =  getPost_downloadForVisualising_v1(iPt)


%iPt = 11;
% 
% addpath('IEEGToolbox');
% addpath('IEEGToolbox/lib');
% addpath('ieeg-cli-1.13.4');
% addpath('ieeg-cli-1.13.4/config');
% addpath('ieeg-cli-1.13.4/lib');

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
LeadTime = 5*60*60;

% Get type 3 seizures = 1
Type3 = 0;

% Time before & after sz in seconds to save
Tbefore = 1*60;
Tafter = 60*60;

% training data time cutoff (days)
start_cutoff = 15*7;
end_cutoff = 60*7;

%% Feature parameters
% FILTERS
NV_filters_EN

%% load information
load(['Portal Annots/' curPt '_Annots']);
load('Portal Annots/portalT0');
trial_t0 = datenum(startDateTime(iPt));

% chron. order
[SzTimes,I] = sort(SzTimes);


SzType = SzType(I);
SzDur = SzDur(I);

% circadian times
SzCirc = trial_t0 + SzTimes/1e6/86400;
SzCirc = datevec(SzCirc);
SzCirc = SzCirc(:,4);
SzDay = SzTimes/1e6/86400;
shift = 1;

%% get seizure index
ISI = diff(SzTimes)/1e6;
ISI = [LeadTime+1 ISI];
if ~Type3
    remove = SzType == 3;
    SzType(remove) = [];
    ISI(remove) = [];
    SzTimes(remove) = [];
    SzCirc(remove) = [];
    SzDur(remove) = [];
    SzDay(remove) = [];
end

% save seizures that have a leading interval of LeadTime and within
% % training period
% %SzDay = ceil(SzTimes/1e6/60/60/24);
% SzInd = find(ISI > LeadTime);

SzDay = ceil(SzTimes/1e6/60/60/24);
training = SzDay > start_cutoff & SzDay < end_cutoff;
SzInd = find(ISI > LeadTime & training);

%% start grabbing data
N = length(SzInd);
%preIctal = zeros(Window,N);
%IctalDropouts = zeros(Window,N);
IctalCirc = zeros(1,N);
segment_length = 5*60;

meanTime = ceil(Fs_actual*segment_length);


timeS = round(Fs_actual*5);

fprintf('%d seizures\n',N)
for n = 1:N
    fprintf('Seizure %d of %d\n',n,N)
    % intialize start time
    t0 = SzTimes(SzInd(n))/1e6 - Tbefore;
    
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
    
    % pre-filter
    Data(isnan(Data(:,1)),:) = 0;
%     Data = filtfilt(filter_wb(1,:),filter_wb(2,:),Data);
%     Data = filtfilt(filter_notch(1,:),filter_notch(2,:),Data);
    
    % need to grab the data in segments
    window = Tbefore + ceil(SzDur(n)) + Tafter;
    
    for nn = 0:shift:(window/shift-shift)
        ind1 = floor(Fs_actual*nn)+1;
        
        try
            curSeg = Data(ind1:ind1+timeS,:);
        catch
            curSeg = Data(ind1:end,:);
        end
        
        if sum(curSeg(:,1).^2) < 1e-16
            % ignore dropout sections
            IctalDropouts(nn/shift+1,n) = 1;
            continue;
        end
        
        
    end % end feature segments
    
    % save the time of day
    IctalCirc(n) = SzCirc(n);
    
%% Displaying example traces
    Fs=400;                     % Approx sample rate of EEG
    NChs = 16;                  % Number of channels
    NSamples = size(Data,1);     % Number of samples
    
    t = 1/Fs:1/Fs:NSamples/Fs;  % Time vector in seconds
    
    %Zero mean the data instead of highpass filter
    ZM_data = Data - repmat(mean(Data,1)',1,NSamples)';
    
    F_Ord = 2;                  % Filter Order
    Fc = 35;                    % Hz, filter cut off, low for visualization
    Wn = Fc/(Fs/2);             % Normalized cutoff
    [b, a] = butter(F_Ord, Wn); % define LP filter
    
    F_ZM_data = filtfilt(b,a,ZM_data);  % Apply lowpass filter
    
    Ch_offset = 400;                                    % Distance between channels
    Offset_vector = Ch_offset:Ch_offset:Ch_offset*NChs; % Offset for each channel
    Offset_Mat = repmat(Offset_vector', 1, NSamples);   % Extend vector to a matrix
    
    OS_F_ZM_data = F_ZM_data + Offset_Mat';             % Add offset to the data
    
    sz_len = SzDur(SzInd(n));                           % Get length of seizure info
    
    sz_start_ind = Tbefore;                             % Determine the start and end of the seizures
    sz_end_ind = (sz_len)+ Tbefore;
    
%% Save variables for later use
    

    save_dir = ['C:/Users/depayne/Desktop/PostIctalExamples/Pt_' num2str(iPt) '/'];    % Directory to save all variables into
    save_path = [save_dir 'Sz_' num2str(n) '.mat'];             % Path to save all variables to
    try
        save(save_path)                                         % Save variables to file
    catch 
        mkdir(save_dir)                                         % Make patient specific folder if not there
        save(save_path)
    
    end

%     csvwrite([parent_path save_path 'Sz_' num2str(n) '.csv'],Data);
    fprintf('%d of %d epochs processed\n',n,N)
    
    pause(10) % give matlab enough time to finish writing the file to disk
    %TESTING!
end % end seizure loop

% csvwrite([parent_path save_path 'SzDropouts.csv'],IctalDropouts);
% csvwrite([parent_path save_path 'SzHour.csv'],IctalCirc);

