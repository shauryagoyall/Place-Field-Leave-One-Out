function place_fields= calculate_place_fields_epochs(x_bins_width,timestamps,varargin)
% INPUTS:
%   x_bin_width: enter width value (2 for fine resolution, 10 for bayesian
%   decoding).
%   timestamps: list of cells with start:stop matrix for each track
% loads list_of_parameters.m, extracted_clusters.mat,extracted_position.mat, extracted_waveform.mat
% uses function  skaggs_information.m

parameters = list_of_parameters;

load('extracted_clusters.mat');
load('extracted_position.mat');
if exist('extracted_waveforms.mat','file')
    load('extracted_waveforms.mat');
else
    disp('no extracted_waveforms.mat file');
    allclusters_waveform=[];
end

if x_bins_width == parameters.x_bins_width_bayesian
    load('extracted_place_fields_BAYESIAN.mat');
elseif x_bins_width == parameters.x_bins_width
    load('extracted_place_fields.mat');
else
    disp('loading default place fields')
    load('extracted_place_fields.mat');
end


%find track positions
for track_id=1:length(position.linear)
    track_indices = ~isnan(position.linear(track_id).linear);
    place_fields.track(track_id).time_window=[min(position.t(track_indices)) max((position.t(track_indices)))];
end

% Run threshold on pyramidal cells: half-width amplitude
if ~isempty(allclusters_waveform)
    PC_indices = [allclusters_waveform.half_width] > parameters.half_width_threshold; % cells that pass treshold of pyramidal cell half width
    pyramidal_cells = [allclusters_waveform(PC_indices).converted_ID];
end

%find mean rate spikes on all tracks / time on all tracks
%(for identifying putative interneurons later in the code)
for track_id = 1:length(position.linear)
    total_time_in_track(track_id) = place_fields.track(track_id).time_window(2)-place_fields.track(track_id).time_window(1);
    for j = 1 : max(clusters.id_conversion(:,1))
        all_spikes(track_id,j) = length(find(clusters.spike_id==j & ...
            clusters.spike_times>place_fields.track(track_id).time_window(1) & ...
            clusters.spike_times<place_fields.track(track_id).time_window(2)));
    end
end
place_fields.mean_rate=sum(all_spikes,1)/sum(total_time_in_track);


%% Place field calculation
time_bin_width=position.t(2)-position.t(1);
for track_id= 1:length(position.linear)
    epochs= timestamps{track_id};
    if ~isempty(epochs)
        % find positions
        pos_idx=[];
        spike_idx= [];
        for this_epoch=1:size(epochs,1)
            pos_idx= [pos_idx find(position.t >= epochs(this_epoch,1) & position.t <= epochs(this_epoch,2))];
            spike_idx= [spike_idx find(clusters.spike_times >= epochs(this_epoch,1) & clusters.spike_times <= epochs(this_epoch,2))'];
        end
        pos= position.linear(track_id).linear(pos_idx);
        spike_times= clusters.spike_times(spike_idx);
        place_fields.track(track_id).time_window= epochs;
        total_time= cumsum(epochs(:,2)- epochs(:,1));
        total_time= total_time(end);
        place_fields.track(track_id).total_time= total_time;

        position_speed = abs(position.v_cm(pos_idx));
        nan_idx= isnan(pos);
        position_speed(nan_idx) = NaN;  %make sure speed is NaN if position is NaN
        position_during_spike = interp1(position.t(pos_idx),pos, spike_times,'nearest'); %interpolates position into spike time
        speed_during_spike = interp1(position.t(pos_idx),position_speed,spike_times,'nearest');

        % Time spent at each x_bin (speed filtered)
        x_bin_edges = 0:x_bins_width:100*position.linear(track_id).length; % forces x_bins to be from 0 to the maximum length (converting from m cm)
        x_bin_centres = [(x_bin_edges(2)-x_bins_width/2):x_bins_width:(x_bin_edges(end-1)+x_bins_width/2)];
        x_hist = time_bin_width.*histcounts(pos(find(position_speed>parameters.speed_threshold_laps...
            & position_speed<parameters.speed_threshold_max)),x_bin_edges);
        place_fields.track(track_id).x_bin_centres = x_bin_centres;
        place_fields.track(track_id).x_bin_edges = x_bin_edges;
        place_fields.track(track_id).x_bins_width = x_bins_width;
        place_fields.track(track_id).dwell_map = x_hist;

        for j=1:max(clusters.id_conversion(:,1))
            % Number of spikes per bin within time window (speed filtered)
            place_fields.track(track_id).spike_hist{j} = histcounts(position_during_spike(find(clusters.spike_id(spike_idx)==j & ...
                speed_during_spike>parameters.speed_threshold_laps &...
                speed_during_spike<parameters.speed_threshold_max)),x_bin_edges);
            place_fields.track(track_id).raw{j} = place_fields.track(track_id).spike_hist{j}./x_hist; % place field calculation
            place_fields.track(track_id).raw{j}(find(isnan(place_fields.track(track_id).raw{j})))=0;  %replace all NaNs with the value 0, to allow smoothing

            % zero bins with 0 dwell time, but make sure no spikes occurred
            non_visited_bins = find(x_hist==0);
            if sum(place_fields.track(track_id).spike_hist{j}(non_visited_bins))>0
                disp('ERROR: x_hist is zero, but spike histogram is not');
            else
                place_fields.track(track_id).raw{j}(non_visited_bins)= 0;
            end
            place_fields.track(track_id).non_visited_bins = non_visited_bins; %NaNs that have been replaced by O

            % Create smoothing filter
            if x_bins_width== parameters.x_bins_width_bayesian
                w= [1 1];  %moving average filter of 2 sample, will be become a filter of [0.25 0.5 0.25] with filtfilt
            else
                w= gausswin(parameters.place_field_smoothing);
            end
            w=w./sum(w); %make sure smoothing filter sums to 1

            % Get place field information
            place_fields.track(track_id).smooth{j}         =    filtfilt(w,1,place_fields.track(track_id).raw{j}); %smooth pl field
            place_fields.track(track_id).centre_of_mass(j) =    sum(place_fields.track(track_id).smooth{j}.*x_bin_centres)/sum(place_fields.track(track_id).smooth{j});  %averaged center
            [place_fields.track(track_id).peak(j) , index] =    max(place_fields.track(track_id).smooth{j}); %peak of smoothed place field and index of peak (center)
            if place_fields.track(track_id).peak(j) ~=0
                if length(index)>1 % very rare exception where you have multiple peaks of same height....
                    index= index(1);
                end
                place_fields.track(track_id).centre(j) = x_bin_centres(index);
                % add delimitation to place field
                firing_thresh = 0.2 * place_fields.track(track_id).peak(j);
                firing_idx = place_fields.track(track_id).smooth{j} > firing_thresh;
                firing_idx_before= firing_idx(1:index);
                firing_idx_after= firing_idx(index:end);
                edge_idx_before = diff(firing_idx_before);
                edge_idx_after = diff(firing_idx_after);
                % find place cell boudaries
                pc_bounds_1 = find(edge_idx_before == 1);
                pc_bounds_2 = find(edge_idx_after == -1)+(index-1);
                if isempty(pc_bounds_1)  % no defined limits
                    pc_bounds_1= 1; % first index
                end
                if isempty(pc_bounds_2)
                    pc_bounds_2= length(place_fields.track(track_id).smooth{j});
                end
                % find nearest threshold crossings to peak (in case threshold is crossed
                % mulitple times)
                place_fields.track(track_id).field_boundaries(j,:)= [max(pc_bounds_1) min(pc_bounds_2)];
            else
                place_fields.track(track_id).centre(j) = NaN;
            end
            place_fields.track(track_id).raw_peak(j)          = max(place_fields.track(track_id).raw{j}); % raw pl field peak
            place_fields.track(track_id).mean_rate_session(j) = length(find(clusters.spike_id==j))/(position.t(end)-position.t(1)); %mean firing rate
            place_fields.track(track_id).mean_rate_track(j)   = sum(place_fields.track(track_id).spike_hist{j})/total_time;
            if place_fields.track(track_id).peak(j) ~=0
                place_fields.track(track_id).half_max_width(j) = x_bins_width*half_max_width(place_fields.track(track_id).smooth{j}); %finds half width of smoothed place field (width from y values closest to 50% of peak)
            else
                place_fields.track(track_id).half_max_width(j) = NaN;
            end
        end

        %calculate skagges information
        %place_fields.track(track_id).skaggs_info= skaggs_information(place_fields.track(track_id));
        
    end

% Sort new cells by place field centre
[~,index] = sort(place_fields.track(track_id).centre);
place_fields.track(track_id).sorted = index;
[~,index] = sort(place_fields.track(track_id).centre(place_fields.track(track_id).good_cells));
place_fields.track(track_id).sorted_good_cells = place_fields.track(track_id).good_cells(index);
[~,index] = sort(place_fields.track(track_id).centre(place_fields.track(track_id).good_cells_LIBERAL));
place_fields.track(track_id).sorted_good_cells_LIBERAL = place_fields.track(track_id).good_cells_LIBERAL(index);

end

if x_bins_width == parameters.x_bins_width_bayesian
    place_fields.good_place_cells= place_fields_BAYESIAN.good_place_cells;
    place_fields.good_place_cells_LIBERAL= place_fields_BAYESIAN.good_place_cells_LIBERAL;
    place_fields.unique_cells= place_fields_BAYESIAN.unique_cells;
    place_fields.unique_cells= place_fields_BAYESIAN.unique_cells;
    place_fields.pyramidal_cells= place_fields_BAYESIAN.pyramidal_cells;
    place_fields.other_cells= place_fields_BAYESIAN.other_cells;
elseif x_bins_width == parameters.x_bins_width
    place_fields.good_place_cells= place_fields.good_place_cells;
    place_fields.good_place_cells_LIBERAL= place_fields.good_place_cells_LIBERAL;
    place_fields.unique_cells= place_fields.unique_cells;
    place_fields.unique_cells= place_fields.unique_cells;
    place_fields.pyramidal_cells= place_fields.pyramidal_cells;
    place_fields.other_cells= place_fields.other_cells;
end






end


function half_width = half_max_width(place_field)
%interpolate place field to get better resolution
new_step_size=0.1;  %decrease value to get finer resolution interpolation of place field
place_field_resampled=interp1(1:length(place_field),place_field,1:new_step_size:length(place_field),'linear');
[peak,index] = max(place_field_resampled); %finds smoothed place field peak firing rate (FR)
for i = index : length(place_field_resampled)
    if place_field_resampled(i)<peak/2 %finds the point after the peak where the FR is half the peak FR
        break;
    end
end
for j = index : -1 : 1 %finds the point before the peak where the FR is half the peak FR
    if place_field_resampled(j)<peak/2
        break;
    end
end
half_width = new_step_size*(i-j); %distance between half-peaks
%(calculated in indicies of original place field, but converted to distance in cm in function above)
end
