clear;
clc;
load('significant_replay_events.mat')
load('extracted_place_fields_BAYESIAN.mat')
load('extracted_position.mat')

% POST_start=max(position.linear(2).timestamps);  % for sleep (actually remote replay)
% POST_end=min(position.linear(3).timestamps);

POST_start=min(position.linear(1).timestamps);  %for awake replay
POST_end=max(position.linear(2).timestamps);


spike_counts=zeros(length(significant_replay_events.track),length(place_fields_BAYESIAN.mean_rate));
event_counts=zeros(length(significant_replay_events.track),length(place_fields_BAYESIAN.mean_rate));

for track_id = 1:length(significant_replay_events.track)
   
    %count cells
    track_events = significant_replay_events.track(track_id).spikes;
    num_events = length(track_events);

    for event_num = 1:num_events
        if significant_replay_events.track(track_id).event_times(event_num)>=POST_start & significant_replay_events.track(track_id).event_times(event_num)<=POST_end
            
            event = track_events(1,event_num);
            for i = 1:length(event{1,1})
                cell = event{1,1}(i,1);
                spike_counts(track_id,cell) = spike_counts(track_id,cell)+1;
            end

            unique_cells = unique(event{1,1}(:,1));
            for i = 1:length(unique_cells)
                cell = unique_cells(i);
                event_counts(track_id,cell) = event_counts(track_id,cell)+1;
            end
        end
    end
    
    % %count events
    % for event_num = 1:num_events
    %     event = track_events(1,event_num);
    %     unique_cells = unique(event{1,1}(:,1));
    %     for i = 1:length(unique_cells)
    %         cell = unique_cells(i);
    %         event_counts(track_id,cell) = event_counts(track_id,cell)+1;
    %     end
    % end

    %cells that did not replay
    fprintf('Good cells on track %d that did not replay ', track_id);
    good_cells = place_fields_BAYESIAN.track(track_id).good_cells;
    num_good_cells = length(good_cells);
    for i = 1:num_good_cells
        cell = good_cells(i);
        if spike_counts(track_id,cell)==0
            disp(cell);
        end
    end

    fprintf("Good cells across any track that did not replay ");
    good_cells = place_fields_BAYESIAN.good_place_cells;
    num_good_cells = length(good_cells);
    for i = 1:num_good_cells
        cell = good_cells(i);
        if spike_counts(track_id,cell)==0
            disp(cell);
        end
    end

    %print most active cells by event
    [sorted_events, sorted_indices] = sort(event_counts(track_id,:), 'descend');
    for i = 1:length(spike_counts(track_id,:))
        cell_number = sorted_indices(i);
        event_count = sorted_events(i);
        if event_count ~= 0
            fprintf('Cell %d: Events %d\n', cell_number, event_count);
        end
    end

    %print most active cells by spike
    cell_id = 1:length(spike_counts);
    [sorted_spikes, sorted_indices] = sort(spike_counts(track_id,:), 'descend');
    for i = 1:length(spike_counts(track_id,:))
        cell_number = sorted_indices(i);
        spike_count = sorted_spikes(i);
        if spike_count ~= 0
            fprintf('Cell %d: spikes %d\n', cell_number, spike_count);
        end
    end

    event_counts_this_track = event_counts(track_id,:);

    replay_counts(track_id).spike_count = sorted_spikes;
    replay_counts(track_id).cell_id = sorted_indices;
    replay_counts(track_id).event_counts = event_counts_this_track(sorted_indices);
end

save('replay_counts_awake.mat', 'replay_counts');
