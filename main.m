%% BIBLIOGRAPHY
% Reinhard Sonnleitner, PhD Thesis - 2017 - Audio Identification via Fingerprinting Achieving Robustness to Severe Signal Modifications 
% Sébastien Fenet, PhD Thesis - 2013 - Audio-Fingerprints And Associated Indexing Strategies For The Purpose Of Large-Scale Audio-Identification
% Avery Li-Chun Wang, Article - 2003 - An Industrial-Strength Audio Search Algorithm

%% CLEAR
clear all
close all
clc
cd('C:/Users/Apeth/Desktop/Projects/Shazam') ;                             %set current folder

%% PARAMETRISATION
param = struct(...
'DSrate' , 4 ,...                                                          %downsampling rate
'filter_order' , 1 ,...                                                    %Butterworth filter order
'frame_size_time' , 64 ,...                                                %window size for spectrogram computing [ms]
'overlap' , 30 ,...                                                        %windows overlaping rate for spectrogram computing [%]
'window_type' , 'hann' ,...                                                %window type (bartlett, blackman, chebwin, hamming, hann, kaiser, triang, rectwin, ...
'c' , 10 ,...                                                              %noise (%max) added to the query sample [%]
'deltat' , 1 ,...                                                          %histogram bin size [s]
'delta_F_tile_hz' , 400 ,...                                               %size of spectrum regions/divisions for maxima searching [Hz]
'delta_T_tile_s' , 0.2 ,...                                                %size of spectrum regions/divisions for maxima searching [s]
'delta_F_max_hz' , 350 ,...                                                %target zone size(frequency) [Hz]
'delta_T_max_s' , 0.3 ,...                                                 %target zone size(time) [s]
'n1' , 13 ,...                                                             %number of bits to encode f(anchor)  
'n2' , 13 ,...                                                             %number of bits to encode f(point)
'n3' , 6 ,...                                                              %number of bits to encode t(point)-t(anchor)
'Tvote' , 50 ,...                                                          %threshold, rate of coherents results needed [%]
'Lsub' , 1 ,...                                                            %length of the sub-queries [s]
'osub' , 50 ,...                                                           %overlap between the successive sub-queries [%]
't_sample' , 5 ,...                                                       %set sample length [s]
'plot', 0 );                                                               %plot audio fingerprint construction

%% DATABASE DIRECTORY
dirName = '.\test_songs\';                                                 %set the directory where songs are

%% BUILD THE DATABASE
[index_dB,SONGS] = build_database(dirName,param);                          %create an index with audio fingerprints from each song

%% RANDOMLY SELECT A SAMPLE
[X_sample,fs_sample,ref] = sample_generation(dirName,param);               %randomly chose a signal from database and add noise

%% DIVIDE THE SAMPLE IN OVERLAPPING FRAMES
[sw,indices] = signalintoframes(X_sample,fs_sample,param) ;                %divide the signal into overlapping frames

%% CREATE THE AUDIO FINGERPRINT of the SAMPLE
index_sample = compute_audiofingerprint_sample(sw,fs_sample,param) ;       %create an index (cell) with the audio fingerprint from the sample

%% INDEX SEARCH
[ID,DT,MAX_HIST] = index_search(index_dB,index_sample,SONGS,param) ;       %for each frame, find the best match in the database index and returns the ID of the song and the difference of time

%% RESULTS ANALYSIS
results_analysis(ID,DT,MAX_HIST,ref,param) ;


%%idea : add density computation for each song and mean
