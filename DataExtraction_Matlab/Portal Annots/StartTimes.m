% set up the start times for the portal

close all
clc
clear

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

DateNum(1) = datenum('10/06/2010, 07:40:34','dd/mm/yyyy, HH:MM:SS');
DateNum(2) = datenum('02/08/2010, 23:59:05','dd/mm/yyyy, HH:MM:SS');
DateNum(3) = datenum('15/11/2010, 02:01:34','dd/mm/yyyy, HH:MM:SS');
DateNum(4) = datenum('10/11/2010, 01:23:15','dd/mm/yyyy, HH:MM:SS'); 
DateNum(5) = datenum('10/05/2011, 02:48:58','dd/mm/yyyy, HH:MM:SS');
DateNum(6) = datenum('08/06/2011, 02:20:27','dd/mm/yyyy, HH:MM:SS');
DateNum(7) = datenum('23/07/2010, 03:17:28','dd/mm/yyyy, HH:MM:SS');
DateNum(8) = datenum('19/11/2010, 02:16:02','dd/mm/yyyy, HH:MM:SS');
DateNum(9) = datenum('27/05/2011, 01:43:00','dd/mm/yyyy, HH:MM:SS');
DateNum(10) = datenum('07/06/2011, 06:49:26','dd/mm/yyyy, HH:MM:SS');
DateNum(11) = datenum('08/07/2010, 06:45:38','dd/mm/yyyy, HH:MM:SS');
DateNum(12) = datenum('08/07/2010, 22:14:18','dd/mm/yyyy, HH:MM:SS');
DateNum(13) = datenum('02/08/2010, 02:47:09','dd/mm/yyyy, HH:MM:SS');
DateNum(14) = datenum('25/11/2010, 00:33:27','dd/mm/yyyy, HH:MM:SS');
DateNum(15) = datenum('06/05/2011, 09:34:10','dd/mm/yyyy, HH:MM:SS');

% DATE TIMES IN UTC
startDateTime = datetime(datevec(DateNum),'TimeZone','UTC');

% MATLAB DATENUMS CONVERTED TO MELBOURNE
startAEST = TimezoneConvert(startDateTime,'UTC','Australia/Melbourne');

save('portalT0','Patient','startDateTime','startAEST');

