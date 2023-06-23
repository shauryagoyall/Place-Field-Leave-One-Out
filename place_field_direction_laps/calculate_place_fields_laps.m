function varargout= calculate_place_fields_laps(varargin)
% options are 'direction':
                    % bidirectional (default)  or unidirectional
% 'grouping':
                % 'each' (default) for individual laps (full or half laps depending on direction option)
                % 'last_laps'
                % cell array (in dev)
% 'size': 
                % size of bins, default is parameters.x_bins_width
% 'n_laps':
                % if last_laps is selected in 'grouping', n_laps = 4 if no number
                % entered (default)
% 'save_option': 0 or 1

parameters= list_of_parameters;
load('lap_times.mat');

default_n_laps=0;

p = inputParser;
addParameter(p,'direction','unidirectional',@ischar);
addParameter(p,'grouping','even_odd',@ischar);
addParameter(p,'size',parameters.x_bins_width,@isnumeric);
addParameter(p,'n_laps',default_n_laps,@isnumeric);
addParameter(p,'save_option',1,@isnumeric);
parse(p,varargin{:});

n_laps= p.Results.n_laps;
if strcmp(p.Results.grouping,'last_laps') && contains('n_laps',p.UsingDefaults)
    n_laps= 4;
end

% build filename from options
str_name= []; str_name2=[];
if strcmp(p.Results.direction,'unidirectional')
    str_name= [str_name '_directional_'];
end
str_name= [str_name 'place_fields'];

if strcmp(p.Results.grouping,'each')
    str_name= [str_name '_all_laps'];
elseif strcmp(p.Results.grouping,'even_odd')
    str_name= [str_name '_even_odd'];
elseif strcmp(p.Results.grouping,'last_laps')
    str_name= [str_name '_comp_last_' num2str(n_laps) '_laps'];
    str_name2= ['place_fields_last_' num2str(n_laps) '_laps'];
end
% if p.Results.size == parameters.x_bins_width_bayesian
%     str_name= [str_name '_BAYESIAN'];
%     str_name2= [str_name2 '_BAYESIAN'];
% elseif p.Results.size == parameters.x_bins_width
%     % nothing
% else % arbitrary size
%     str_name= [str_name num2str(p.Results.size) 'cm'];
%     str_name2= [str_name2 num2str(p.Results.size) 'cm'];
% end
% str_name= [str_name '.mat'];
% str_name2= [str_name2 '.mat'];

    for this_track=1:length(lap_times)
        typ= mod(length(lap_times(this_track).lap),2);
        if strcmp(p.Results.direction,'bidirectional')
            if typ % odd number of laps
                start_idx= 1:2:length(lap_times(this_track).lap)-1;
                stop_idx= 2:2:length(lap_times(this_track).lap)-1;
            else % even number of laps
                start_idx= 1:2:length(lap_times(this_track).lap);
                stop_idx= 2:2:length(lap_times(this_track).lap);
            end
        elseif strcmp(p.Results.direction,'unidirectional')
                start_idx= 1:1:length(lap_times(this_track).lap);
                stop_idx= 1:1:length(lap_times(this_track).lap);
        else
            disp('option not recognised');
            keyboard;
        end
        if strcmp(p.Results.grouping,'last_laps')
            % if n_laps
            last_laps{this_track}= [lap_times(this_track).start(start_idx(end-(n_laps-1):end))' lap_times(this_track).end(stop_idx(end-(n_laps-1):end))'];
            previous_laps{this_track}= [lap_times(this_track).start(start_idx(1:end-n_laps))' lap_times(this_track).end(stop_idx(1:end-n_laps))']; 
        elseif strcmp(p.Results.grouping,'each')
            previous_laps{this_track}= [lap_times(this_track).start(start_idx(1:end-n_laps))' lap_times(this_track).end(stop_idx(1:end-n_laps))']; 
        elseif strcmp(p.Results.grouping,'even_odd')
            even_laps{this_track}= [lap_times(this_track).start(start_idx(2:2:end))' lap_times(this_track).end(stop_idx(2:2:end))'];
            odd_laps{this_track}= [lap_times(this_track).start(start_idx(1:2:end))' lap_times(this_track).end(stop_idx(1:2:end))']; 
        end
    end

    % calculate fields for last n laps
    if n_laps>0 
        place_fields_last_laps= calculate_place_fields_epochs(p.Results.size,last_laps);
        % plot_place_fields(place_fields_last_laps);
        save(str_name2,'place_fields_last_laps','-v7.3');
    end

    % pad arrays with zeros
    if ~strcmp(p.Results.grouping,'even_odd')
        length_laps= cell2mat(cellfun(@(x) size(x,1),previous_laps,'UniformOutput',false));
        max_laps= max(cell2mat(cellfun(@(x) length(x),previous_laps,'UniformOutput',false)));
        to_pad= max_laps - length_laps;
        previous_laps= cellfun(@(x,y) padarray(x,[y 0],'post'),previous_laps,num2cell(to_pad),'UniformOutput',false);

        for this_track=1:length(lap_times)
        %     disp(['Track...' num2str(this_track)])
            for this_lap= 1:length_laps(this_track)
        %         disp(['        lap...' num2str(this_lap)])
                lap_ts= cell(1,length(lap_times));
                lap_ts{this_track}= previous_laps{this_track}(this_lap,:);
                fld_temp= calculate_place_fields_epochs(p.Results.size,lap_ts);

                % calculate place fields for each lap
                place_field_laps(this_track).(['lap' num2str(this_lap)])= fld_temp;
                place_field_laps(this_track).(['lap' num2str(this_lap)]).track= fld_temp.track(this_track);

            end
        end
    elseif strcmp(p.Results.grouping,'even_odd')
        place_fields_even= calculate_place_fields_epochs(p.Results.size,even_laps);
        place_fields_odd= calculate_place_fields_epochs(p.Results.size,odd_laps);
    end
    % 
    % % save under different names depending on options
    % if p.Results.save_option
    %     save("extracted_place_field_directional",'place_field_laps','-v7.3');
    %     disp("saved")
    % end
if p.Results.save_option
         save("extracted_place_field_laps_even",'place_fields_even','-v7.3');
         save("extracted_place_field_laps_odd",'place_fields_odd','-v7.3');
         disp("saved")
end    


    if ~strcmp(p.Results.grouping,'even_odd')
        varargout{1}= place_field_laps;
        if n_laps>0
            varargout{2}= place_fields_last_laps;
        end
    elseif strcmp(p.Results.grouping,'even_odd')
        varargout{1}= place_fields_even;
        varargout{2}= place_fields_odd;
    end
      
end