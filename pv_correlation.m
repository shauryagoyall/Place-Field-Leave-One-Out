clear;
clc;
close all;

load("extracted_place_field_laps");
load("extracted_place_fields");

%%% add the cell number from clusters here 
good_cells = place_fields.good_place_cells;

population_data = place_field_laps;
for track_id = 1:length(place_field_laps)
    fieldNames = fieldnames(place_field_laps(track_id));
    for i = 1:numel(fieldNames)
       lap_id = fieldNames{i};  % Get the current field name
       current_lap_place_field = place_field_laps(track_id).(lap_id); 

       population_stack = [];
       if ~isempty(current_lap_place_field)
           for j = 1:length(good_cells)
               population_stack = vertcat(population_stack,current_lap_place_field.track.smooth{1,good_cells(j)});
           end
       else
           break;
       end
       population_data(track_id).(lap_id)=population_stack;
    end

end

correlated_values = place_field_laps;
for track_id = 3:length(place_field_laps)
    fieldNames = fieldnames(place_field_laps(track_id));
    for i = 1:numel(fieldNames)-1
       lap_id = fieldNames{i} ;  % Get the current field name
       current_lap_population = population_data(track_id).(fieldNames{i});
       next_lap_population = population_data(track_id).(fieldNames{i+1});

       lap_correlation = zeros(1, size(current_lap_population,2) );

       if ~isempty(next_lap_population)
           for j = 1:size(current_lap_population,2)
               corr_matrix = corrcoef(current_lap_population(:,j),next_lap_population(:,j));
               lap_correlation(1,j) = min( corr_matrix(1, 2), 1);
           end
       else
           correlated_values(track_id).(lap_id) = [];
           break;
       end
         correlated_values(track_id).(lap_id) = [];
         correlated_values(track_id).(lap_id).mean = mean(lap_correlation);
         correlated_values(track_id).(lap_id).distribution = lap_correlation;
    end

    figure;
    correlations = [];
    for i = 1:numel(fieldNames)-1
        if ~isempty(correlated_values(track_id).(fieldNames{i}))
            correlations = [correlations correlated_values(track_id).(fieldNames{i}).mean];
        else
            break;
        end
    end
    plot(correlations)
    xlabel("Laps")
    ylabel("Correlation")
    xlim([0, 16])
    ylim([0 1])
    title("Track "+num2str(track_id))

end

