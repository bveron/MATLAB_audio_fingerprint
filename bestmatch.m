function [id,dt,MAX_hist] = bestmatch(index_dB,index_sample_frame,SONGS,param)
    T = [] ;
    for  ii = 1 : size(index_sample_frame,1)
        
        a1 = index_dB(:,1)==index_sample_frame(ii,1) ;
        a2 = index_dB(:,2)==index_sample_frame(ii,2) ;
        a3 = index_dB(:,3)==index_sample_frame(ii,3) ;
        matches = find( a1.*a2.*a3 == 1) ;
        
        for jj = 1 : size(matches,1)
%         	if index_dB(jj,1:3) == index_sample_frame(ii,1:3)
%         matches = find(ismember(index_dB(:,1:3),index_sample_frame(i,1:3)),1) ;
            T = [T ; index_dB(matches(jj),5) index_dB(matches(jj),4)-index_sample_frame(ii,4)] ;
%             end
            
        end
    end
    
    DTmax = [] ;
    Amplitude = [] ;

    for ii = 1 : size(SONGS,1)
%         figure ;
        tancc = T(:,1)==ii ;
%         subplot(1,size(SONGS,1),ii)
%         hold on
%         histogram(T(tancc,2),floor(SONGS(ii,2)/param.deltat))
        HIST = histogram(T(tancc,2),floor(SONGS(ii,2)/param.deltat));
        ind_max = find(HIST.Values==max(HIST.Values)) ;
        DTmax = [DTmax ; HIST.BinEdges(ind_max(1))'] ;
        Amplitude =  [Amplitude ; HIST.Values(ind_max(1))'] ;
    end
    
    id = find(Amplitude==max(Amplitude)) ;
    id = id(1) ;
    dt = DTmax(id) ;
    MAX_hist = max(Amplitude) ;

%     match_max_time = find(HIST.Values == max(HIST.Values));
%     time_song = HIST.BinEdges(match_max_time(1));
% 
%     MATCH_MAX_TIME(i) = match_max_time(1) ;
%     TIME_SONG(i) = time_song;
%     T = cell(size(index_sample_frame,1),2) ;
%     tic
% 	for i = 1 : size(index_sample_frame,1) 
%         matches = find(index_dB(:,1) == index_sample_frame(i,1)) ;
%         T{i,1} = i ;
%         T{i,2} = matches ;
%     end
%     
%     T_length = cellfun('length',T) ;
%     T_match = T_length(:,1) .* T_length(:,2) ;
%     all_matches = sum(T_match) ;
%     %there is probably a better way to do this, T_match could be useful
%     tanchor = zeros(all_matches,1) ;
%     id = zeros(all_matches,1) ;
%     coun = 0 ;
%     for i = 1 : size(T,1)
%         match_indDB = T{i,2};
%         if length(match_indDB)~=0
%             for j = 1 : length(match_indDB)
%                 coun = coun + 1 ;
%                 tanchor(coun) = index_dB(match_indDB(j),2);
%                 id(coun) = index_dB(match_indDB(j),3) ;
% %                 R(coun) = match_indDB(j) ;
% %                 H(coun) = [T{i,1}] ;
%             end
%         end
%     end
%     
% %     MATCHES = [id,dt] ;
%     
%     MATCH_MAX_TIME = zeros(size(SONGS,1),1);
%     TIME_SONG = zeros(size(SONGS,1),1);
%     %%recherche du plus grand nombre de matches par seconde de chanson
%     for i = 1 : size(SONGS,1)
% %         figure ; histogram(dt(id==i),floor(SONGS(i,2)/param.deltat))
%         HIST = histogram(tanchor(id==i),floor(SONGS(i,2)/param.deltat)) ;
%         match_max_time = find(HIST.Values == max(HIST.Values));
%         time_song = HIST.BinEdges(match_max_time(1));
%         
%         MATCH_MAX_TIME(i) = match_max_time(1) ;
%         TIME_SONG(i) = time_song;
%         
%     end
%     
%     id = find(MATCH_MAX_TIME == max(MATCH_MAX_TIME)) ;
%     id = id(1) ;
%     
% %     disp('---------------------------------------------------------------------------------------')
% %     disp(['Best match is song number: ',num2str(indMAX)])
%     
% 	tanch = TIME_SONG(id) ;
    
    
    
    
%     nb_song = mode(id) ;
    
%     figure ; histogram(DT,round(Ltot)) ;
%     xlabel('Time diffence betwen query and maximum matches')
%     title('Time difference between query and DB')
%     
%     XX = histogram(Time,round(Ltot)) ;
%     match_max_time = find(XX.Values == max(XX.Values));
%     time_song = XX.BinEdges(match_max_time);
%     disp([num2str(time_song),' seconds from the beginning'])
      
    
    
%     position_matches_index_dB = R ;
    
%     figure; hist(position_matches_index_dB,round(Ltot/param.deltat));
%     xlabel('Second of database')
%     ylabel('Match ammount')
%     title('Matched couple indexdB position')
    
%     figure; hist(index_dB(position_matches_index_dB,4),round(Ltot/param.deltat));
%     xlabel('Second of database')
%     ylabel('Matches ammount')
%     title('Matched anchor time position')
    
%     figure; hist(index_dB(position_matches_index_dB,3),index_dB(end,3));  
%     xlabel('Song ID')
%     ylabel('Matches ammount')
%     title('Song ID')
    
%     match_max_id = find(hist(index_dB(position_matches_index_dB,3),index_dB(end,3)) == max(hist(index_dB(position_matches_index_dB,3),index_dB(end,3))));
    
%     ind_max_in_index_dB = match_max_id*size(index_dB,1)/round(Ltot/param.deltat) ;
    
%     id_song = index_dB(round(ind_max_in_index_dB),5) ;
%     disp('---------------------------------------------------------------------------------------')
%     disp(['Best match is song number: ',num2str(match_max_id)])
%     count = 0 ;
%         
%     for i = 1 : size(T,1)
%         if size(T{i,2},1)~=0
%             A = T{i,2};
%             for j = 1 : length(A)
% %                 if index_dB(A(j),3) == match_max_id
%                    count = count + 1 ;
%                    Time(count) = index_dB(A(j),2) - index_dB(T{i,1},2);
% %                 end
%             end
%         end
%     end
%     figure ; histogram(Time,round(Ltot)) ;
%     xlabel('Time diffence betwen query and maximum matches')
%     title('Time difference between query and DB')
%     
%     XX = histogram(Time,round(Ltot)) ;
%     match_max_time = find(XX.Values == max(XX.Values));
%     time_song = XX.BinEdges(match_max_time);
%     disp([num2str(time_song),' seconds from the beginning'])
    
    
    
    
end
