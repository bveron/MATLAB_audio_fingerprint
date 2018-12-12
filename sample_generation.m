function [X_sample,fs_sample,ref] = sample_generation(dirName,param) 
    
    files = dir('C:\Users\Apeth\Desktop\Projects\Shazam\test_songs\**\*.mp3');                  
    files = {files.name}';
    
    tira1 = round(((rand(1,1)*(length(files)-1))+1)) ;
    disp('---------------------------------------------------------------------------------------')
    disp(['Song chosen is number: ',num2str(tira1)]) ;
    filetoload = char(files(tira1));
    [x,fs_sample] = audioread([dirName,filetoload]) ;
    Xtemp = mean(x,2) ;
    
    L_song_ind = length(Xtemp) ;
    L_song_s = L_song_ind / fs_sample ;
    
    ind_max = floor((L_song_s-param.t_sample)*fs_sample) ;
    
    tira2 = round(((rand(1,1)*(ind_max-1))+1)) ;
    disp([num2str(tira2/fs_sample),' seconds from the beginning']) ;
    indXs = tira2 ;
    t_sample_ind = floor(param.t_sample*fs_sample) ;
    
    X_sample = Xtemp(indXs:indXs+t_sample_ind-1);
    
    ma = max(X_sample) ;
    noise = (param.c/100)*ma*2*(rand(length(X_sample),1)-1);
%     SNR = snr(X_sample,noise); %signal-to-noise ration in dB
    X_sample = X_sample + noise ;
    
    clear x Xtemp
    
    ref.id = tira1 ;
    ref.time = tira2/fs_sample ;
      
end
