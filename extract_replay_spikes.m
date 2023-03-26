%function out=extract_replay_spikes

clear;
clc;
close all;

import decoding.calculate_estimated_position

out(1).decoded_replay=[];
out(2).decoded_replay=[];
load('significant_replay_events.mat')
load('extracted_place_fields_BAYESIAN.mat')
load('extracted_position.mat')

POST_start=max(position.linear(2).timestamps);  %I think this will extract the start and end time of sleep POST 1....
POST_end=min(position.linear(3).timestamps);


good_place_cells=place_fields_BAYESIAN.good_place_cells;

for track_id=1:2

    for j=1:length(place_fields_BAYESIAN.track(track_id).raw) %number of place cells on track
        place_fields{j}=place_fields_BAYESIAN.track(track_id).raw(j);
    end
    for j=1:length(significant_replay_events.track(track_id).spikes) %number of replay events
      
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
%                     %position_bins=linspace(0,1,length(good_place_cells));
%                   
%                     [position,estimated_position_time,estimated_position_interp, estimated_position]=calculate_estimated_position(t,t0,bin_width,place_fields(place_fields_BAYESIAN.track(1).good_cells),event_spike_times,event_spike_id,position_bins,good_place_cells);
%                     %place_fields(:,good_place_cells)
%                     %%%%%%%%%%FIGURE 2: Bayesian decoding
%                     figure;
%                     % %imagesc(estimated_position_time,position_bins,position)
%                     imagesc(estimated_position_time,position_bins,position)
%                     a=colormap(bone);
%                     colormap(flipud(a));
%                     axis xy
%                     hold on
%                     % plot(t,estimated_position_interp_actual,'r','LineWidth',2)
%                     plot(t0,estimated_position,'r','LineWidth',2)
%                     %plot(t,estimated_position_interp,'r','LineWidth',2)
%                     title(['Bayesian decoding- Neuron' num2str(left_spike_id) ' left out- only good'])
%                     ylabel('Position (normalized)')
%                     xlabel('Time(s)')
%                     
%                     colorbar
                    
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
                    

                    figure;
                    subplot(1 ,2 ,1);
                    x_ind = [];
                    for i = 1:length(event_spike_id)
                        if ismember(event_spike_id(i),place_fields_BAYESIAN.track(1).good_cells)
                            x_ind = cat(2,x_ind,place_fields{1,event_spike_id(i)}{1,1}');
                        end
                    end
                    imagesc(flipud(x_ind));
                    colormap(hot);
                    title("Track 1 good only")
                    
                    subplot(1 ,2 ,2);
                    selected = place_fields(event_spike_id);
                    mat = [];
                    for i = 1:length(selected)
                        mat = cat(2,mat,selected{1,i}{1,1}');
                    end   
                    imagesc(flipud(mat));
                    colormap(hot);
                    %ylabel([0:2:20])
                     title("All in event")
                   
                    %%plot the neurons that are skipped here to get a sanity check
                    %left_out_spike_times = cell_spike_times;
                    %left_out_spike_position = position_bins(left_spike_id);
                    %plot(left_out_spike_times,left_out_spike_position,'r.','linewidth',1000);                      
                end
            end
        end
    end
end


