function index_sample = compute_audiofingerprint_sample(sw,fs_sample,param) 
    
    index_sample = cell(1,size(sw,1)) ;
    file_id = 0 ;

    for ii = 1 : size(sw,1)
        clear index_sample_frame
        [index_sample_frame] = compute_audiofingeprint(sw(ii,:),fs_sample,param,file_id) ;
        
        index_sample_frame(:,4) = index_sample_frame(:,4) + (ii-1)*param.Lsub*(1-param.osub/100) ;

        index_sample{ii} = index_sample_frame ;
        
    end
    
end
