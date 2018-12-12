function ref_found = results_analysis(ID,DT,MAX_HIST,ref,param)

    cou = 0 ;

    ri = ID ;
    si = DT ;
    marge = (param.deltat/2) ;

    si_focus_min = abs(round(si - marge)) ;
    si_focus_max = abs(round(si + marge)) ;


    UN = unique(abs(round(si))) ;

    MATCH = sum(diff(si)<1) ;%%margin : 1 second

    
    figure;
    subplot(3,1,1)
    plot(ID,'b+');
    title('best ID match /frame')
    hold on;
    subplot(3,1,2)
    plot(DT,'b+');
    hold on;
    refline(0,ref.time)
    title('best DT match /frame')
    hold on;
    subplot(3,1,3)
    plot(MAX_HIST,'b+')
    title('MAX HIST')
    
    if MATCH > length(ID)*(param.Tvote/100)
        disp('MATCH FOUND')
        kz = find(MAX_HIST == max(MAX_HIST)) ;
        kz = kz(1) ;
        disp(['ID: ',num2str(ID(kz))])
        disp(['time: ',num2str(DT(kz))])
        
        hold on;
        subplot(3,1,1)
        plot(kz,ID(kz),'r+');
        subplot(3,1,2)
        plot(kz,DT(kz),'r+');
        subplot(3,1,3)
        plot(kz,MAX_HIST(kz),'r+');
        
        ref_found.id = ID(kz) ;
        ref_found.time = DT(kz) ;
        
    else
        
        disp('MATCH NOT FOUND')
    
        ref_found.id = 0 ;
        ref_found.time = 0 ;
        
    end
    
end