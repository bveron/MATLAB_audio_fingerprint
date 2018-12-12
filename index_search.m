function [ID,DT,MAX_HIST] = index_search(index_dB,index_sample,SONGS,param) 
   
    Lsize = length(index_sample) ;
    ID = zeros(Lsize,1) ;
    DT = zeros(Lsize,1) ;
    MAX_HIST = zeros(Lsize,1) ;
    
    for ii = 1 : Lsize
        index_sample_frame = index_sample{ii} ;
        [id,dt,MAX_hist] = bestmatch(index_dB,index_sample_frame,SONGS,param) ;
        ID(ii) = id ;
        DT(ii) = dt ;
        MAX_HIST(ii) = MAX_hist ;
    end

    
end
