function plot_place_fields(varargin)
%plots place fields on each track, sorted by each of the tracks (2 tracks=4
%plots, 3 tracks=9 plots). first plot is a heat map of place field, second
%plot is a firing rate vs. position plot
%
%leave input empty if you want to use previous extracted place fields saved
%as a .mat file



if isempty(varargin)
    place_fields = importdata('extracted_place_field_laps_odd.mat'); 
else
    place_fields= cell2mat(varargin);
end

% plot_color_line= arrayfun(@(x) {repmat([0.5 0.5 0.5],length(place_fields.track(x).sorted_good_cells),1)},1:length(place_fields.track));
% plot_color_line{1}(place_fields.track(1).sorted_good_cells== 97,1)= 1; plot_color_line{1}(place_fields.track(1).sorted_good_cells== 97,2:3)= 0;
% plot_color_line{2}(place_fields.track(2).sorted_good_cells== 97,1)= 1; plot_color_line{2}(place_fields.track(2).sorted_good_cells== 97,2:3)= 0;
% plot_color_line{1}(place_fields.track(1).sorted_good_cells== 29,3)= 1; plot_color_line{1}(place_fields.track(1).sorted_good_cells== 29,1:2)= 0;
% plot_color_line{2}(place_fields.track(2).sorted_good_cells== 29,3)= 1; plot_color_line{2}(place_fields.track(2).sorted_good_cells== 29,1:2)= 0;

plot_color_line =[139,0,0]/255;

% %%% PLOT WITH IMAGESC
% c=1;
% figure
% y_vector=[];
% for i=3:4
%     for j=1:2
%         matrix=[];
%         normalized_matrix=[];
%         subplot(length(place_fields.track),length(place_fields.track),c)
%         for k=1:length(place_fields.track(j).sorted_good_cells)
%             matrix(k,:)=place_fields.track(i).smooth{place_fields.track(j).sorted_good_cells(k)};
%             normalized_matrix(k,:)=(matrix(k,:)-min(matrix(k,:)))/(max(matrix(k,:))-min(matrix(k,:)));
%         end
%         y_vector= [y_vector, 1.5*i-1];
%         imagesc(normalized_matrix);
%         colormap(jet)
%         hold on
%         ylabel('cell id');
%         yt=place_fields.track(j).sorted_good_cells;
%         set(gca,'yticklabel',yt);
%         xlabel('linearized position')
%         title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(i)]}]);
%         axis xy
%         c=c+1;
%     end
% end

%%% PLOT WITH CURVES
figure;
c=1;
for kk=3:4
    for j=1:2
        y_vector=[];
        color_index=1+floor(place_fields.track(kk).peak*2);
        color_index(find(color_index>64))=64;
        field_color=colormap(jet);
        for ii=1:length(place_fields.track(j).sorted_good_cells)
            %plot sorted
            matrix=[];
            normalized_matrix=[];
            matrix(ii,:)=place_fields.track(kk).smooth{place_fields.track(j).sorted_good_cells(ii)};
            normalized_matrix(ii,:)=(matrix(ii,:)-min(matrix(ii,:)))/(max(matrix(ii,:))-min(matrix(ii,:)));
            subplot(length(place_fields.track)-2,length(place_fields.track)-2,c)
            plfield_row= normalized_matrix(ii,:)+(1.5*ii-1);
            plot(1:length(plfield_row),plfield_row,'k'); hold on;
            x2 = [1:length(plfield_row), fliplr(1:length(plfield_row))];
            inBetween = [(1.5*ii-1)*ones(size(plfield_row)), fliplr(plfield_row)];
            fill(x2, inBetween,field_color(color_index(ii),:)); 
            %fill(x2, inBetween,[139,0,0]/255);
            y_vector= [y_vector, 1.5*ii-1];
        end
        xlim([0 size(normalized_matrix,2)+2]);
        ylim([0 max(y_vector)+1.2]);
        yt=place_fields.track(j).sorted_good_cells;
        set(gca,'ytick',y_vector);
        set(gca,'yticklabel',yt);
        ylabel('Unit ID');
        xlabel('sorted linearized position (bins)');
        c=c+1;
        title([{['place cells on track ' num2str(j)]} ; {['sorted by track ' num2str(kk)]}]);
    end
end

end