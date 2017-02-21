function [] = getInterictal_EN(iPt)

%iPt = 11;

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

temp = ['Inter_' num2str(iPt)];
parent_path = 'C:/Users/depayne/Desktop/NVData/Pt13_09_01_17/';
mkdir(parent_path,temp);
save_path = [temp '/'];

%% link to portal session
curPt = Patient{iPt};
patient = IEEGSession(['NVC1001_' curPt '_2'],login,pword);

%% parameters
Fs_actual = patient.data.sampleRate;
Fs = 400;  % for filtering
iCh = 1:16;

% Get type 3 seizures = 1
Type3 = 0;

% Time in seconds to save
DataWin = 30*60;
% gap between seizure for interictal
InterWin = 5*60*60;
Tbefore = 30*60;

% training data time cutoff (days)
start_cutoff = 2*7;
end_cutoff = inf*7;

shift = 1;                          % shift by 1 s at a time
timeS = round(Fs_actual*5);

% FILTERS
NV_filters_EN

%% load data
% load([curPt '_DataInfo']);
% trial_t0 = datenum(MasterTimeKey(1,:));
load(['Portal Annots/' curPt '_Annots']);
load('Portal Annots/portalT0');
trial_t0 = datenum(startDateTime(iPt));

% chron. order
[SzTimes,I] = sort(SzTimes);
SzType = SzType(I);
SzDur = SzDur(I);
SzInd = SzInd(I);

% circadian times
SzCirc = trial_t0 + SzTimes/1e6/86400;
SzCirc = datevec(SzCirc);

SzDay = SzTimes/1e6/86400;
SzCirc = SzCirc(:,4);

%% get seizure index
ISI = diff(SzTimes)/1e6;
ISI = [0 ISI];
if ~Type3
    remove = SzType == 3;
    SzType(remove) = [];
    ISI(remove) = [];
    SzTimes(remove) = [];
    SzCirc(remove) = [];
    SzDay(remove) = [];
end

% find places where we have at least InterWin hours free of seizure before AND
% after
train_inter = find(ISI > 2*InterWin);
if length(train_inter) == 0
    return
end

train_inter(train_inter == 1) = [];

% this is where we can start interictal windows from (InterWin hours after the seizure before the one with the long ISI)
InterInd = ceil(SzTimes(train_inter-1)/1e6) + InterWin;   % in s
InterCirc = mod(SzCirc(train_inter-1) + InterWin/60/60,24);
% length the interictal window can extend for (InterWin hours before next seizure)
InterLength = floor((ISI(train_inter)-2*InterWin));  % in s

% get rid of times outside training period
training = InterInd > (start_cutoff*24*60*60) & InterInd <= (end_cutoff*24*60*60);
InterLength = InterLength(training);
InterInd = InterInd(training);
InterCirc = InterCirc(training);

InterIctalTimes = [];
InterCircTimes = [];
for m = 1:length(InterInd)
    Time = InterInd(m):DataWin:DataWin*floor((InterInd(m)+InterLength(m))/DataWin);
    CircTime = InterCirc(m):InterCirc(m)+length(Time)-1;    % HACK ONLY WORKS IF DATAWIN IS ONE HOUR
    CircTime = mod(CircTime,24);
    InterIctalTimes = [InterIctalTimes Time];
    InterCircTimes = [InterCircTimes CircTime];
end

%% get rid of times that aren't in the troughs
% load([save_path 'SzProb']);
% lowSeizure = find(SzProb < prctile(SzProb,50))-1;
% valid = ismember(InterCircTimes,lowSeizure);
% InterIctalTimes = InterIctalTimes(valid);
% InterCircTimes = InterCircTimes(valid);
% InterIctalTimes = [InterIctalTimes InterIctalTimes + Tbefore];  % hack to get the extra half an hour windows
% InterCircTimes = [InterCircTimes InterCircTimes];

%%  get the same amount of interictal data
load(['C:\Users\depayne\Desktop\NVData\Pt13_09_01_17\Sz_' num2str(iPt) '\SzHour.csv'])
N = length(SzHour);
randomInd = randperm(length(InterIctalTimes));
segment_length = Tbefore;

nseg = 0;
interIctal = zeros(Tbefore,N);
interIctalCirc = zeros(1,N);
interIctalDropouts = zeros(Tbefore,N);
tic

fprintf('%d seizures\n',N)
n=0;
while nseg < N
    n = n + 1;
    fprintf('Seizure %d of %d\n',n,N)
    if isempty(randomInd)
        display('no more data')
        break;
    end
    
    t0 = InterIctalTimes(randomInd(1));  % interictal window time
    circT = InterCircTimes(randomInd(1));
    randomInd(1) = [];
    
    % get the data from the portal
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
    
    % check dropouts
    if sum(isnan(Data(:,1))) > 0.25*length(Data)
        display('too many dropouts');
        continue;
    end
    Data(isnan(Data(:,1)),:) = 0;
    
    % we can use the segment
    nseg = nseg+1;
    % pre-filter
    Data = filtfilt(filter_wb(1,:),filter_wb(2,:),Data);
    Data = filtfilt(filter_notch(1,:),filter_notch(2,:),Data);
     % need to grab the data in segments
    for nn = 0:shift:(Tbefore/shift-shift);
        ind1 = floor(Fs_actual*nn)+1;
        
        try
            curSeg = Data(ind1:ind1+timeS,:);
        catch
            curSeg = Data(ind1:end,:);
        end
        
        if sum(curSeg(:,1).^2) < 1e-16
            % ignore dropout sections
            interIctalDropouts(nn/shift+1,nseg) = 1;
            continue;
        end
        
        
    end % end feature segments
    
    % save the time of day
    interIctalCirc(nseg) = SzCirc(nseg);
    
    csvwrite([parent_path save_path 'Inter_' num2str(nseg) '.csv'],Data);
    fprintf('%d of %d epochs processed\n',nseg,N)
end % end seizure loop
disp('all done')
csvwrite([parent_path save_path 'InterDropouts.csv'],interIctalDropouts);
csvwrite([parent_path save_path 'InterHour.csv'],interIctalCirc);save