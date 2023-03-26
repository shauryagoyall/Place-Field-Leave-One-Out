%Bayesian decoding 
%function [position,estimated_position_time,estimated_position_interp, estimated_position]=calculate_estimated_position(t,t0,bin_width,place_fields_BAYESIAN,spike_times,spike_id,position_bins,good_place_cells)
function decoded=calculate_estimated_position(t,t0,bin_width,place_fields_BAYESIAN,spike_times,spike_id,position_bins,good_place_cells)
    %disp(length(t0))    
  for track_id = 1:2
        
        %%%% this can be optimized for faster%%%%%%%%%%5
        for j=1:length(place_fields_BAYESIAN.track(track_id).raw) %number of place cells on track
            decoded(track_id).place_fields{j}=place_fields_BAYESIAN.track(track_id).raw(j);
        end
      place_fields(track_id,:)=decoded(track_id).place_fields(:,good_place_cells);   
  end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
  for j=1:length(t0)
            for track_id = 1:2
                temp_position(track_id,:)=reconstruct(t0(j),bin_width,place_fields(track_id,:),spike_times,spike_id,position_bins,good_place_cells)';
            end
            temp_position = renormalize(temp_position, position_bins);
 
            for track_id = 1:2
                decoded(track_id).position(:,j)=temp_position(track_id,:);
            end
  end
  
  for track_id = 1:2      
        for j=1:length(t0)
            [probability(j),position_index(j)]=max(decoded(track_id).position(:,j));
        end
  
        decoded(track_id).estimated_position=position_bins(position_index);
        index=find(probability>0.2);
        decoded(track_id).estimated_position_time=t0+bin_width/2;
        %decoded(track_id).estimated_position_interp=interp1(decoded(track_id).estimated_position_time(index),decoded(track_id).estimated_position(index),t,'linear');
        decoded(track_id).estimated_position_interp=[];
  end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function position=reconstruct(t0,bin_width,place_field,spike_times,spike_id,position_bins,good_place_cells)
        global parameters
        t1=t0+bin_width;
        total_spike_count=0;
        product_of_place_fields=ones(1,length(position_bins));
        sum_of_place_fields=zeros(1,length(position_bins));
        number_of_cells=length(place_field);
        parameters.bayesian_threshold=10.^(log2(number_of_cells)-log2(400)); % small value multiplied to all values to get rid of zeros

        for k=1:number_of_cells

            n=length(find(spike_times>=(t0) & spike_times<=(t1)  & spike_id==good_place_cells(k))); %count # of spikes in bin
            %disp(n)
            
            place_field{1,k}{1,1}(find(place_field{1,k}{1,1}<parameters.bayesian_threshold))=parameters.bayesian_threshold;  %avoid zeros
            single_place_field=interp1(position_bins,place_field{1,k}{1,1},position_bins,'linear');
     
            
            if size(single_place_field,1)>size(single_place_field,2)
                single_place_field=single_place_field';
            end

            product_of_place_fields=product_of_place_fields.*((single_place_field).^(n'*ones(1,length(position_bins)))); %% what is this
            sum_of_place_fields=sum_of_place_fields+single_place_field;
            total_spike_count=total_spike_count+n;

        end

        position=product_of_place_fields.*(exp(-bin_width*sum_of_place_fields)); %bin width is one therfore denom is 1  ;
        %position=position./(sum(position,2)*ones(1,length(position_bins)));
 end
 
 function position=renormalize(position,position_bins)
        for track_id = 1:2
            position(track_id,:)=position(track_id,:)./((sum(position(1,:),2)+sum(position(2,:),2))*ones(1,length(position_bins)));
        end
 end
 