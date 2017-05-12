% This script visualises data downloaded using
% getPost_downloadForVisualising_v1.m


function [] = getPost_Visualise(PatientN, SeizureN)


    load(['C:/Users/depayne/Desktop/PostIctalExamples/Pt_' num2str(PatientN) '/Sz_' num2str(SeizureN) '.mat']);
    
    GuiFigure(OS_F_ZM_data, t, 400, Ch_offset, sz_start_ind*Fs, sz_end_ind*Fs)



end