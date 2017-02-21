% SET FILTERS FOR NV FEATURES

Ford = 2;                                   % filter order

% Wideband Filter
Wc1 = 1;
Wc2 = 100;
W1 = Wc1/(Fs/2);
W2 = Wc2/(Fs/2);
[b,a] = butter(Ford,[W1 W2],'bandpass');
filter_wb = [b;a];

% Notch Filter
Wc1 = 45;
Wc2 = 55;
W1 = Wc1/(Fs/2);
W2 = Wc2/(Fs/2);
[b,a] = butter(Ford,[W1 W2],'stop');
filter_notch = [b;a];