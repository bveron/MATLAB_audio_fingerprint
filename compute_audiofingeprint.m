function index = compute_audiofingeprint(X,fsx,param,file_id)

    %% X PARAMETERS DEFINITION
     LX = length(X) ;
     dtX = 1/fsx ;                                                          %sampling time
     NfftX = 2.^nextpow2(length(X)) ;
     dfX = 1/(NfftX*dtX) ;
     fX = (-NfftX/2 : (NfftX/2)-1)*dfX ;
     tX = (0:1:LX-1)*dtX ;
     
%     FX = ((NfftX/2)-1)*dfX;
%     TX = LX*dtX ;
%     Lref = length(X)*dtX ;

    %% NORMALIZATION
    X = (X-mean(X)) ./ max(X) ;                                            %mean removal and maximum division

    %% FILTERING before downsampling
    Wn = (floor(fsx/2)/param.DSrate)/fsx;                                  %normalized frequency, depends on the downsampling rate for Shannon-Nyquist condition
    [b,a] = butter(param.filter_order,Wn) ;                                %calculation of Butterworth filters coefficients (alternative filters will be implemented)
    Xfil = filter(b,a,X) ;                                                 %apply low-pass filtering
        
    %% DOWNSAMPLING
    Xfil_DS = Xfil(1:param.DSrate:end) ;                                   %downsampling with DSrate
    tXfil_DS = tX(1:param.DSrate:end) ;
    
    if param.plot == 1  
        figure;
        subplot(2,2,1);
        plot(tXfil_DS,X(1:param.DSrate:end),'b')
        hold on
        plot(tXfil_DS,Xfil_DS,'r')
        legend('signal', 'filtered signal')
        title('Decimated signal before and after filtering')
        xlabel('Time (s)')
        ylabel('Amplitude')
        xlim([0 max(tXfil_DS)])
        pause(0.0001)
        
        subplot(2,2,2);
        plot(fX,20*log10(abs(fftshift(fft(X(1:param.DSrate:end),NfftX)))),'b')
        hold on
        plot(fX,20*log10(abs(fftshift(fft(Xfil_DS,NfftX)))),'r');
        legend('signal spectrum', 'filtered signal spectrum')
        title('Spectrum before and after filtering')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (dB)')
        xlim([0 max(fX)])
        pause(0.0001)
    end
    
    %% XDS PARAMETERS DEFINITION
    fsx_DS = fsx / param.DSrate ;                                          %sampling frequency of the downsampled song
%     L_DS = length(Xfil_DS) ;
    dt_DS = 1/fsx_DS ;                                                     %sampling time of the downsampled song
%     Nfft_DS = 2.^nextpow2(length(Xfil_DS)) ;
%     df_DS = 1/(Nfft_DS*dt_DS) ;
%     f_DS = (-Nfft_DS/2 : (Nfft_DS/2)-1)*df_DS ;
%     t_DS = (0:1:Nfft_DS-1)*dt_DS ;
%     F_DS = ((Nfft_DS/2)-1)*df_DS;
%     T_DS = L_DS*dt_DS ;


    %% WINDOWING
    frame_size_samples = floor((param.frame_size_time*1e-3)/dt_DS) ;       %setting frame size, in sample
    noverlap = floor((param.overlap/100)*frame_size_samples) ;             %setting overlap size, in sample
    
    WDW = window(param.window_type,frame_size_samples)' ;                  %window calculation, depends on window type and frame size

    %% XW PARAMETERS DEFINITION
    fsx_W = fsx / param.DSrate ;                                           %sampling frequency of the windowed signal
%     L_W = frame_size_samples ;
    dt_W = 1/fsx_W ;                                                       %sampling time of the windowed signal
    Nfft_W = 2.^nextpow2(frame_size_samples) ;                             %number of samples for FFT calculation in spectrogram computing
    df_W = 1/(Nfft_W*dt_W) ;                                               %frequency step of the FFT of the windowed signal
    f_W = (-Nfft_W/2 : (Nfft_W/2)-1)*df_W ;                                %frequency array of the FFT of the windowed signal
%     t_W = (0:1:Nfft_W-1)*dt_W ;
%     F_W = ((Nfft_W/2)-1)*df_W;
%     T_W = L_W*dt_W ;
%     amount_frames = (length(Xfil_DS)-frame_size_samples)/(frame_size_samples-noverlap) ;

    %% SPECTROGRAM CALCULATION
    [SpecX,Fx,Tx] = spectrogram(Xfil_DS,WDW,noverlap, ...
        f_W((Nfft_W/2)+1:end),fsx_W) ;                                     %spectrogam calculation, for Xfil_DS signal, with WDW window, noverlap overlap, for positive frequencies, with fsx_W sampling frequency
    SpecX = abs(SpecX) ;                                                   %getting the spectrum modulus

    %% XSPEC PARAMETERS DEFINITION
    dT_specX = Tx(2)-Tx(1) ;                                               %spectrogram time step
    dF_specX = Fx(2)-Fx(1) ;                                               %spectrogram frequency step
%     F_specx = Fx(end) ;

    %% SEARCH FOR MAXIMA
    delta_F_tile_ind = floor(param.delta_F_tile_hz / dF_specX) ;           %vertical (frequency) size of a tile in the spectrogram, in sample
    delta_T_tile_ind = floor(param.delta_T_tile_s / dT_specX) ;            %horitonzal (time) size of a tile in the spectrogram, in sample
%     SpecX_bin = SpecX ;

    am_rect_t = floor(size(SpecX,2)/delta_T_tile_ind) ;                    %number of horizontal tiles fitting in the spectrogram
    am_rect_f = floor(size(SpecX,1)/delta_F_tile_ind) ;                    %number of vertical tiles fitting in the spectrogram

%     tile = zeros(delta_F_tile_ind,delta_T_tile_ind) ;

    amount_max = am_rect_t*am_rect_f ;                                     %total number of tiles fitting in the spectrogram, a.k.a. number of maximums

    row = zeros(amount_max,1) ;                                            %initialize the array with maxima vertical (frequency) indexes
    col = zeros(amount_max,1) ;                                            %initialize the array with maxima horizontal (time) indexes

    count = 0 ;
    
    if param.plot == 1    
        subplot(2,2,3); 
        imagesc(Tx,Fx,log(SpecX) ); %plot the log spectrum
        colormap(gray);
        set(gca,'YDir', 'normal');
        hold on
        plot(col*dT_specX,row*dF_specX,'r+')
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        pause(0.0001)
    end
    
    %% SPECTROGRAM DIVIDING
    for i = 1 : am_rect_t
        for j = 1 : am_rect_f

            tile = SpecX( 1+(j-1)*delta_F_tile_ind : j*delta_F_tile_ind , 1+(i-1)*delta_T_tile_ind : i*delta_T_tile_ind ) ;
            [kmaxrecF,kmaxrefT] = find(tile==max(max(tile))) ;
    %         find(SpecX(1+(j-1)*delta_F_tile_ind:j*delta_F_tile_ind,1+(i-1)*delta_T_tile_ind:i*delta_T_tile_ind)==max(SpecX(1+(j-1)*delta_F_tile_ind:j*delta_F_tile_ind,1+(i-1)*delta_T_tile_ind:i*delta_T_tile_ind))) ;
            if size(kmaxrecF,1) > 1
                kmaxrecF = kmaxrecF(1) ;
                kmaxrefT = kmaxrefT(1) ;
            end
               
            count = count + 1 ;
            row(count) = kmaxrecF + (j-1)*delta_F_tile_ind ; 
            col(count) = kmaxrefT + (i-1)*delta_T_tile_ind ;
% 
            if param.plot == 1 
                hold on
                rect = rectangle('Position',[(1+(i-1)*delta_T_tile_ind)*dT_specX (1+(j-1)*delta_F_tile_ind)*dF_specX param.delta_T_tile_s param.delta_F_tile_hz]); 

                hold on
                plot(col(count)*dT_specX,row(count)*dF_specX,'r+')
                pause(0.001)
            end
                        
        end
    end

    IND = sub2ind(size(SpecX),row,col) ;
    IND_sorted = sort(IND) ;
    [row,col] = ind2sub(size(SpecX),IND_sorted) ;

%     Jx = col ;
%     Ix = row ;

    %% TARGET ZONE DEFINITION
    
    amount_maximum = length(row) ;
%     count2 = 0 ; 

    %%TARGET ZONE SEARCHING AND INDEXING

    delta_F_max_ind = floor(param.delta_F_max_hz / dF_specX) ;
    delta_T_max_ind = floor(param.delta_T_max_s / dT_specX) ;        
    
    count3 = 0 ;
    
    if param.plot == 1 
        subplot(2,2,4); 
        imagesc(Tx,Fx,log(SpecX) ); %plot the log spectrum
        colormap(gray);
        set(gca,'YDir', 'normal');
        hold on
        plot(col*dT_specX,row*dF_specX,'r+')
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
    end
        
    for ii = 1 : amount_maximum
        kTZf = find((abs(row(ii+1:end)-row(ii)))<delta_F_max_ind)+ii ;
        kTZt = find((abs(col(ii+1:end)-col(ii)))<delta_T_max_ind)+ii ;
        kTZ = intersect(kTZf,kTZt,'stable') ;
        n_matches = length(kTZ) ;
        count3 = count3 + n_matches ;
    end
     
    index = zeros(count3,5) ;
    count4 = 0 ;
    for i = 1 : amount_maximum
    
    if param.plot == 1 
        hold on ;
        plot(col*dT_specX,row*dF_specX,'r+')
        hold on ;
        plot(col(i)*dT_specX,row(i)*dF_specX,'b+')
    end
    
        kTZf = find((abs(row(i+1:end)-row(i)))<delta_F_max_ind)+i ;
        kTZt = find((abs(col(i+1:end)-col(i)))<delta_T_max_ind)+i ;

        kTZ = intersect(kTZf,kTZt,'stable') ;
        
    if param.plot == 1         
        hold on
        plot(col(kTZf)*dT_specX,row(kTZf)*dF_specX,'g+')
        hold on
        plot(col(kTZt)*dT_specX,row(kTZt)*dF_specX,'y+')
        hold on
        plot(col(kTZ)*dT_specX,row(kTZ)*dF_specX,'k+')
        pause(0.001)
    end
        
        fmax = Fx(end) ;
        deltaTmax = param.delta_T_max_s ;
        
        if isempty(kTZ) == 0
            points_in_target = [row(kTZ) col(kTZ)];
            n_matches = length(kTZ) ;
            for jj = 1 : n_matches
                count4 = count4 + 1 ;
                f1 = row(i)*dF_specX ;
                f2 = points_in_target(jj,1)*dF_specX ;
                deltat = (points_in_target(jj,2)-col(i))*dT_specX ;
%                 f1t = floor((f1/fmax)*((2^(param.n1))-1)) ;
%                 f2t = floor((f2/fmax)*((2^(param.n2))-1)) ;
%                 deltatt = floor((deltat/deltaTmax)*((2^(param.n3))-1)) ;
%                 if param.delta_T_max_s/((2^param.n3)-1) > 2*dT_specX
%                     disp('Issue while encoding deltaTtilde')
%                 end
%                 ktilde = f1t*(2^(param.n2+param.n3))+f2t*(2^param.n3)+deltatt ;
%                 ktilde = str2double(dec2bin(ktilde)) ;
                ind4 = col(i)*dT_specX ;
                index(count4,:) = [deltat,f1,f2,ind4,file_id] ;
            end
        end
    end 
   
    
end
