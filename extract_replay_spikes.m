%function out=extract_replay_spikes

clear;
clc;
%close all;

import decoding.calculate_estimated_position

load('significant_replay_events.mat')
load('extracted_place_fields_BAYESIAN.mat')
load('extracted_position.mat')

POST_start=max(position.linear(2).timestamps);  % for sleep (actually remote replay)
POST_end=min(position.linear(3).timestamps);

% POST_start=min(position.linear(1).timestamps);  %for awake replay
% POST_end=max(position.linear(2).timestamps);

good_place_cells=place_fields_BAYESIAN.good_place_cells;


for track_id=1:2
    for temp =1:length(place_fields_BAYESIAN.track(track_id).raw)
        out(track_id,temp).decoded_replay=[];
        %out(track_id,temp).spike_times=[];
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
                    left_spike_times=cell_spike_times;
                                        
                    t = significant_replay_events.time_bin_edges; %time array for analysis
                    bin_width=.02; % 20ms bin size for replay
                    %the start and end time of replay
                    t0 = spike_times(1)-.05:bin_width:spike_times(end)+.05;
                    position_bins=place_fields_BAYESIAN.track(track_id).x_bin_centres; %20 from 0 to 200
                    decoded=calculate_estimated_position(t0,bin_width,place_fields_BAYESIAN,event_spike_times,event_spike_id,position_bins,good_place_cells,left_spike_times);
                    out(track_id,left_spike_id).decoded_replay=[out(track_id, left_spike_id).decoded_replay decoded(track_id)];
                    %out(track_id,left_spike_id).spike_times=[out(track_id,left_spike_id).spike_times left_spike_times'];
                    
%            %%%%% TO PLOT THE DECODED CELLS         
%                     for track_num=1:2
%                         place_fields=decoded(track_num).place_fields;
%                         figure;
% 
%                         % %imagesc(estimated_position_time,position_bins,position)
%                         imagesc(decoded(track_num).estimated_position_time,position_bins,decoded(track_num).position)
%                         a=colormap(bone);
%                         colormap(flipud(a));
%                         axis xy
%                         hold on
%                         % plot(t,estimated_position_interp_actual,'r','LineWidth',2)
%                         plot(t0,decoded(track_num).estimated_position,'r','LineWidth',2)
%                         %plot(t,estimated_position_interp,'r','LineWidth',2)
%                         title(['Neuron' num2str(left_spike_id) ' left out'])
%                         ylabel('Position (normalized)')
%                         xlabel('Time(s)')
%                         colorbar
%                     
%                         figure;
%                         subplot(1 ,2 ,1);
%                         x_ind = [];
%                         for i = 1:length(event_spike_id)
%                             if ismember(event_spike_id(i),place_fields_BAYESIAN.track(track_num).good_cells)
%                                 x_ind = cat(2,x_ind,place_fields{1,event_spike_id(i)}{1,1}');
%                             end
%                         end
%                         imagesc(flipud(x_ind));
%                         colormap(hot);
%                         title("Track 1 good only")
% 
%                         subplot(1 ,2 ,2);
%                         selected = place_fields(event_spike_id);
%                         mat = [];
%                         
%                         for i = 1:length(selected)
%                             mat = cat(2,mat,selected{1,i}{1,1}');
%                         end   
%                         
%                         imagesc(flipud(mat));
%                         colormap(hot);
%                         title("All in event")
%   

%%%% Write code to plot the neuron that was skipped on top of the above
%%%% plot to get a sanity check of when the neuron fires. It can be plotted
%%%% as a straight line parallel to y axis at the time points when it fires
%%%% as we are currently assuming that we do not know the place field of
%%%% the neuron

                end
            end
        end
    end
end

save('decoded_replay_sleep','out') %% change sleep to awake for awake replay
