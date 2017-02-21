clear
close all
clc

iPts = [11];
% iPts = 14;
for n = 1:length(iPts)
    getSeizures_EN(iPts(n))
%     getInterictal_EN(iPts(n))
end
