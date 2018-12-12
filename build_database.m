function [index_dB,SONGS] = build_database(dirName,param)

    list_of_songs = dir('.\test_songs\**\*.mp3');                        %get all .mp3 files names form dirName directory        
    list_of_songs = {list_of_songs.name}';
    index_dB = [] ;                                                        %initialize the index of the database (concatenation of fingerprints from each song)
    SONGS = zeros(length(list_of_songs),2) ;                               %initialize a matrix catching songs Id and length, in seconds, for histogram calculation in bestmatch.m
    
    for file_id = 1 : length(list_of_songs)                                %for each song in the directory
        
        clc
        disp('---------------------------------------------------------------------------------------')
        disp(['Audio-fingerprinting the database: ',num2str(file_id),' / ',num2str(length(list_of_songs))])
        filetoload = char(list_of_songs(file_id));                         
        [x,fsx] = audioread([dirName,filetoload]) ;
        X = mean(x,2) ;                                                    %stereo to mono conversion
        clear x                                                            %quickly cleaning memory space
        SONGS(file_id,1) = file_id ;
        SONGS(file_id,2) = length(X)*(1/fsx) ;       
        
        index = compute_audiofingeprint(X,fsx,param,file_id) ;             %compute the audiofingerprint from the current song
                
        index_dB = [index_dB;index] ;                                      %add the index from the current song to the database index
        
    end  
end
        
        
        
    




