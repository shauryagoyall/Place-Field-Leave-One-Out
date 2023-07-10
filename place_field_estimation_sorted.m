clear;
clc;
load('decoded_replay_sleep.mat') %%change sleep to awake for awake replay
load('extracted_place_fields_BAYESIAN.mat')
load('replay_counts.mat')

for track_id = 1:2
    relevant_decoded = out(track_id,:);
    good_place_cells=replay_counts(track_id).cell_id  ;
    relevant_decoded = relevant_decoded(good_place_cells);
    position_bins=place_fields_BAYESIAN.track(track_id).x_bin_centres;
    never_replay_cell = [];
    for i = 70:length(good_place_cells)
        decoded_replay = relevant_decoded(i).decoded_replay;
        cell_id = good_place_cells(i);
        if isempty(decoded_replay)
            never_replay_cell=[never_replay_cell cell_id];
            continue
        end
        spikes_in_position_bin = zeros(1,length(position_bins));
        %spikes_in_position_bin_test = zeros(1,length(position_bins));
        time_in_position_bin = zeros(1,length(position_bins));
        for j = 1:length(decoded_replay)
            left_spike_times = decoded_replay(j).spike_times;
            estimated_position_time = decoded_replay(j).estimated_position_time;
            
            %% spikes in position bin
            for k=1:length(left_spike_times) %for each spike time, find the interval it is between and add 1 to the number of spikes in that interval bin
                for l=1:length(estimated_position_time)
                    if left_spike_times(k) >= estimated_position_time(l)-0.01 & left_spike_times(k) <estimated_position_time(l)+0.01
                       if decoded_replay(j).include_bin(l) ~= 0  
                            ind = find(position_bins==decoded_replay(j).estimated_position(l));
                            spikes_in_position_bin(ind) = spikes_in_position_bin(ind) + 1 ;
                       end
                       
%%%%%%%%%%% to check if the no spike time bins are being skipped                       
%                        ind = find(position_bins==decoded_replay(j).estimated_position(l));
%                        spikes_in_position_bin_test(ind) = spikes_in_position_bin_test(ind) + 1 ;

                    end
                end
            end
        

            %% time in position bin
            for k = 1:length(decoded_replay(j).estimated_position)
                position_index = find(position_bins==decoded_replay(j).estimated_position(k));
                time_in_position_bin(position_index) = time_in_position_bin(position_index) + 1;
            end
            
        end
        
%%%%%%%%%%% to check if the no spike time bins are being skipped         
%         if ~isequal(spikes_in_position_bin,spikes_in_position_bin_test)
%             fprintf("track id %.2f cell id  %.2f \n", track_id, cell_id);
%         end
        
        %% Place Field
        for t=1:length(position_bins)
            frequency(t) = spikes_in_position_bin(t)/time_in_position_bin(t); %10 as time interval is 0.1s
        end
          
%         %% Plots
        figure;
        subplot(3,1,1);
        bar(position_bins,spikes_in_position_bin)
        xticks(position_bins)
        xlabel("Estimated Position")
        ylabel("Number of Spikes")
        title("Number of spikes in a position bin")

        subplot(3,1,2);
        bar(position_bins,time_in_position_bin)
        xticks(position_bins)
        xlabel("Estimated Position")
        ylabel("Time spent")
        title("Time in a position bin")
        
        subplot(3,1,3);
        hold on;
        plot(position_bins,rescale(frequency,0,1));
        actual_place = rescale(decoded_replay(1).place_fields{1, cell_id}{1, 1},0,1 );
        plot(position_bins, actual_place ,'r' ); 
        xticks(position_bins);
        xlabel("Estimated Position");
        ylabel("Firing Frequency");
        title("Estimated Place field "+num2str(cell_id));
        legend("Estimated","Actual");
        hold off;
         
    end%%% put a flag here and run in debug mode to see each plot one by one
    disp("The cells that never replayed are"+ num2str(never_replay_cell))
end

