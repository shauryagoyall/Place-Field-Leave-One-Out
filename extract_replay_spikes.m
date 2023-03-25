%function out=extract_replay_spikes

clear;
clc;
close all;


out(1).decoded_replay=[];
out(2).decoded_replay=[];
load('significant_replay_events.mat')
load('extracted_place_fields_BAYESIAN.mat')
load('extracted_position.mat')

POST_start=max(position.linear(2).timestamps);  %I think this will extract the start and end time of sleep POST 1....
POST_end=min(position.linear(3).timestamps);

%count = 0;
good_place_cells=place_fields_BAYESIAN.good_place_cells;
%disp(length(good_place_cells))
for track_id=1:2
    %disp(count)
    %count = 0;
    for j=1:length(place_fields_BAYESIAN.track(track_id).raw) %number of place cells on track
        place_fields{j}=place_fields_BAYESIAN.track(track_id).raw(j);
    end
    for j=75:length(significant_replay_events.track(track_id).spikes) %number of replay events
      
        if significant_replay_events.track(track_id).event_times(j)>POST_start & significant_replay_events.track(track_id).event_times(j)<POST_end   %only analyze replay during POST sleep

            spike_id=significant_replay_events.track(track_id).spikes{j}(:,1);
            spike_times=significant_replay_events.track(track_id).spikes{j}(:,2);
            
            for k=1:length(good_place_cells)
                cell_id=good_place_cells(k);
                
                if ~isempty(find(spike_id==cell_id)) %if not empty where spike id=cell id means if current good cell spikes in replay event then perfrom bayesian decoding
                    event_spike_id=spike_id;
                    event_spike_times=spike_times;
                    index=find(spike_id==cell_id);
                    cell_spike_times=spike_times(index);   %times cell fires in replay event
                    event_spike_id(index)=[];  %remove spike times for the cell analyzed
                    event_spike_times(index)=[];%remove spike times for the cell analyzed
                    left_spike_id = cell_id; %the cell id that was ignored
                    
                    t = significant_replay_events.time_bin_edges; %time array for analysis
                    bin_width=.02; % 20ms bin size for replay
                    %t0=t(1):bin_width:(t(end)-bin_width); %the start and end time of replay
                    t0 = significant_replay_events.track(track_id).event_times(j)-.1:bin_width:significant_replay_events.track(track_id).event_times(j)+.1;
                    
                    position_bins=place_fields_BAYESIAN.track(track_id).x_bin_centres; %20 from 0 to 200
                    %position_bins=linspace(0,1,length(good_place_cells));
                  
                    [position,estimated_position_time,estimated_position_interp, estimated_position]=calculate_estimated_position(t,t0,bin_width,place_fields(place_fields_BAYESIAN.track(1).good_cells),event_spike_times,event_spike_id,position_bins,good_place_cells);
                    %place_fields(:,good_place_cells)
                    %%%%%%%%%%FIGURE 2: Bayesian decoding
                    figure;
                    % %imagesc(estimated_position_time,position_bins,position)
                    imagesc(estimated_position_time,position_bins,position)
                    a=colormap(bone);
                    colormap(flipud(a));
                    axis xy
                    hold on
                    % plot(t,estimated_position_interp_actual,'r','LineWidth',2)
                    plot(t0,estimated_position,'r','LineWidth',2)
                    %plot(t,estimated_position_interp,'r','LineWidth',2)
                    title(['Bayesian decoding- Neuron' num2str(left_spike_id) ' left out- only good'])
                    ylabel('Position (normalized)')
                    xlabel('Time(s)')
                    
                    colorbar
                    
                    [position,estimated_position_time,estimated_position_interp, estimated_position]=calculate_estimated_position(t,t0,bin_width,place_fields(:,good_place_cells),event_spike_times,event_spike_id,position_bins,good_place_cells);
                    %
                    %%%%%%%%%%FIGURE 2: Bayesian decoding
                    figure;
                    % %imagesc(estimated_position_time,position_bins,position)
                    imagesc(estimated_position_time,position_bins,position)
                    a=colormap(bone);
                    colormap(flipud(a));
                    axis xy
                    hold on
                    % plot(t,estimated_position_interp_actual,'r','LineWidth',2)
                    plot(t0,estimated_position,'r','LineWidth',2)
                    %plot(t,estimated_position_interp,'r','LineWidth',2)
                    title(['Bayesian decoding- Neuron' num2str(left_spike_id) ' left out all neruon'])
                    ylabel('Position (normalized)')
                    xlabel('Time(s)')
                    
                    colorbar
                    
                    %scatter plot x is time y is position 
                    %take the cell id from event spike id
                    %take the location from place_fields_BAYESIAN.track(1).sorted_good_cells
                    % take the time from event spike times
                    %plot
                    %scatter(x,y)
                    figure;
                    x_ind = [];
                    for i = 1:length(event_spike_id)
                        if ismember(event_spike_id(i),place_fields_BAYESIAN.track(1).good_cells)
                            x_ind = cat(2,x_ind,place_fields{1,event_spike_id(i)}{1,1}');
                        end
                    end
                    imagesc(flipud(x_ind));
                    colormap(hot);
                    title("trck 1 good only")
                    
                    
                    figure;
                    selected = place_fields(event_spike_id);
                    mat = [];
                    for i = 1:length(selected)
                        mat = cat(2,mat,selected{1,i}{1,1}');
                    end   
              
                    imagesc(flipud(mat));
                    colormap(hot);
                     title("all in event")
                   
                    %xlim([7.82 7.83])
                    %spike_index = find(spike_id==s)

                    %%plot the neurons that are skipped here to get a sanity check
                    %left_out_spike_times = cell_spike_times;
                    %left_out_spike_position = position_bins(left_spike_id);
                    %plot(left_out_spike_times,left_out_spike_position,'r.','linewidth',1000);
                       
                    % 2) determine which position bin has highest posterior probabilityfor each time bin
                    % 3) calculate number of spikes in each time bin
                    % 4) calculate place field using all cumulative time
                    % bins
                
                %peak probailities below thrieshold (10% of overall prob (about P>=.1)
                %skip time bins if no peak posiition 
                end
            end
        end
    end
end


%Bayesian decoding 
function [position,estimated_position_time,estimated_position_interp, estimated_position]=calculate_estimated_position(t,t0,bin_width,place_field,spike_times,spike_id,position_bins,good_place_cells)
    %disp(length(t0))    
    
        for j=1:length(t0)
            position(:,j)=reconstruct(t0(j),bin_width,place_field,spike_times,spike_id,position_bins,good_place_cells)';
            %disp(size(position))
        end
        
        for j=1:length(t0)
            [probability(j),position_index(j)]=max(position(:,j));
        end
        estimated_position=position_bins(position_index);
        index=find(probability>0.2);
        estimated_position_time=t0+bin_width/2;

        estimated_position_interp=interp1(estimated_position_time(index),estimated_position(index),t,'linear');
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
          
            %what position bins do i take initially ? place fields are 20
            %so i assume 0 to 200 track divided by 20 and then the place
            %field rate corresponds to the rate at the center
            
            %place field not place field bayesian
            %need smoothening of placefield buefroe putting in 
            
            %too few spikes to decode
      
            
            if size(single_place_field,1)>size(single_place_field,2)
                single_place_field=single_place_field';
            end

            product_of_place_fields=product_of_place_fields.*((single_place_field).^(n'*ones(1,length(position_bins)))); %% what is this
            sum_of_place_fields=sum_of_place_fields+single_place_field;
            total_spike_count=total_spike_count+n;

        end

        
        position=product_of_place_fields.*(exp(-bin_width*sum_of_place_fields))/bin_width; %bin width is one therfore denom is 1  ;
        position=position./(sum(position,2)*ones(1,length(position_bins)));

    end