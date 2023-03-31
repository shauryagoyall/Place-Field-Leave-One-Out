clear;
clc;
load('decoded_replay_new.mat')
load('extracted_place_fields_BAYESIAN.mat')
good_place_cells=place_fields_BAYESIAN.good_place_cells;


for track_id = 1:2
    relevant_decoded = out(track_id,:);
    relevant_decoded = relevant_decoded(good_place_cells);
    position_bins=place_fields_BAYESIAN.track(track_id).x_bin_centres;
    for i = 1:length(good_place_cells)
        decoded_replay = relevant_decoded(i).decoded_replay;
        cell_id = good_place_cells(i);
        
        spikes_in_position_bin = zeros(1,length(position_bins));
        time_in_position_bin = zeros(1,length(position_bins));
        for j = 1:length(decoded_replay)
            left_spike_times = decoded_replay(j).spike_times;
            estimated_position_time = decoded_replay(j).estimated_position_time;
            
            %% spikes in position bin
            for k=1:length(left_spike_times) %for each spike time, find the interval it is between and add 1 to the number of spikes in that interval bin
                for l=1:length(estimated_position_time)
                    %if left_spike_times(k) >= estimated_position_time(l)-0.01 & left_spike_times(k) <estimated_position_time(l)+0.01
                    if left_spike_times(k) >= estimated_position_time(l) & left_spike_times(k) <estimated_position_time(l+1)
                        ind = find(position_bins==decoded_replay(j).estimated_position(l));
                        spikes_in_position_bin(ind) = spikes_in_position_bin(ind) + 1 ;
                    end
                end
            end
        
        
%         %linkaxes(ax,'x');
%         spikes_in_position_bin = zeros(1, length(position_bins));
%         spike_index = find(spikes_in_time_bin ~=0);
%         for i=1:length(spike_index)
%             position_index = find( estimated_position(spike_index(i)) == position_bins     );
%             spikes_in_position_bin(position_index) = spikes_in_position_bin(position_index) + 1;
%         end

            %% time in position bin
            for k = 1:length(decoded_replay(j).estimated_position)
                position_index = find(position_bins==decoded_replay(j).estimated_position(k));
                time_in_position_bin(position_index) = time_in_position_bin(position_index) + 1;
            end
        end
        %% Spikes in a position bin
        figure;
        subplot(3,1,1);
        bar(position_bins,spikes_in_position_bin)
        xticks(position_bins)
        xlabel("Estimated Position")
        ylabel("Number of Spikes")
        title("Number of spikes in a position bin")

        %% Time in a position bin
        subplot(3,1,2);
        bar(position_bins,time_in_position_bin)
        xticks(position_bins)
        xlabel("Estimated Position")
        ylabel("Time spent")
        title("Time in a position bin")

        %% Place Field
        for i=1:length(position_bins)
            frequency(i) = 10*spikes_in_position_bin(i)/time_in_position_bin(i); %10 as time interval is 0.1s
        end
        subplot(3,1,3);
        hold on;
        plot(position_bins,frequency);
        %if max(decoded_replay(1).place_fields{1, cell_id}{1, 1}) ~= 0
        actual_place = rescale(decoded_replay(1).place_fields{1, cell_id}{1, 1},0,max(frequency) );
        plot(position_bins, actual_place ,'r' ); 
       % end
        xticks(position_bins);
        xlabel("Estimated Position");
        ylabel("Firing Frequency");
        title("Estimated Place field");
        hold off;
        
    end
    
%% Calculating the place field again
end

