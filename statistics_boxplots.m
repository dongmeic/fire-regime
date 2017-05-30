% close all; clear; clc;

tic

%%% statistics_boxplot.m : Box plot
%% Input : several maps in the same resolution
%% Output: figure and table
%% Description:
%% Created by:

basedir = '/Volumes/dongmeic/fire/output/revision/results/';
% basedir = 'variables/';

%% File names.
filenames = {
	     %'mean_annual_number_of_fires.tif'; 
         	    'mean_annual_area_burned.tif'; 
                      %'maximum_fire_size.tif';
          'mean_annual_active_fire.tif';
          %'mean_annual_number_of_fires_cv.tif';
          'mean_annual_area_burned_cv.tif'; 
	     'mean_annual_active_fire_cv.tif'; 
	     'fire_season_duration.tif'; 
         %'length_of_fire_period.tif';
         %'length_of_fire_period_BA.tif';
             %'peak_start.tif'; 
             'peak_month.tif';
             'fire_radiative_power.tif'; 
              'Gini_index.tif';
	     'percentage_forest_affected.tif';
	     'percentage_savanna_affected.tif';
	     'percentage_grassland_affected.tif';
	     'percentage_cropland_affected.tif';
        %'percentage_forest_affected_BA.tif';
       %'percentage_savanna_affected_BA.tif';
       %'percentage_grassland_affected_BA.tif';
       %'percentage_cropland_affected_BA.tif';       
}

file_name = {
	     %'mean annual number of fires';
         'mean annual burned area';
                %'maximum fire size';
          'mean annual active fire';
         %'CV of annual number of fires';
         'CV of annual burned area';
	     'CV of annual active fire';
	     'fire season duration';
         %'length of fire period';
		     %'peak start';
         'peak month';
         'fire radiative power';
         	     'Gini index';
	     'percentage forest affected';
	     'percentage savanna affected';
	     'percentage grassland affected';
	     'percentage cropland affected';	     
}

%labels = {'MANF'; 'MAFD'; 'CVNF'; 'CVFD'; 'FSD'; 'FPM'; 'FRP'; 'GI'; 'PFA'; 'PSA'; 'PGA'; 'PCA'};
labels = {'MAAB'; 'MAFD'; 'CVAB'; 'CVFD'; 'FSD'; 'FPM'; 'FRP'; 'GI'; 'PFA'; 'PSA'; 'PGA'; 'PCA'};
[map ref] = geotiffread([basedir filenames{1}]);
[m n] = size(map);
len = length(filenames);

for k = 1:len
  map = geotiffread([basedir filenames{k}]);
  if ~all(size(map) == [m n])
    fprintf('[Error]: size is mismatched!\n');
  end
  InputData(:, k) = map(:);
end

%% Boxplot
figure('Units','normalized')
for k = 1:len
  map = geotiffread([basedir filenames{k}]);
  map2=double(map);
  subplot(3,4,k,'replace');
  boxplot(map2(map2 > 0),'labels',[num2str(k) '-' labels{k}],'OutlierSize',4);
end
print([basedir 'clusters/boxplot_variables_' date '.png'],'-dpng','-r500')


% Table
if 1
  filename = ([basedir 'clusters/Statistics_variables_' date '.csv']);
  x = zeros(10,18);
  for k = 1:len
    map = geotiffread([basedir filenames{k}]);
    map = double(map);
    x(k, 1) = k;
    x(k, 2) = prctile(map(map>0),2.5);
    x(k, 3) = prctile(map(map>0),5);
    x(k, 4) = prctile(map(map>0),10);
    x(k, 5) = prctile(map(map>0),25);
    x(k, 6) = prctile(map(map>0),30);
    x(k, 7) = prctile(map(map>0),35); 
    x(k, 8) = prctile(map(map>0),50);
    x(k, 9) = prctile(map(map>0),65);
    x(k, 10) = prctile(map(map>0),70);
    x(k, 11) = prctile(map(map>0),75);
    x(k, 12) = prctile(map(map>0),90);
    x(k, 13) = prctile(map(map>0),95);
    x(k, 14) = prctile(map(map>0),97.5);
    x(k, 15) = min(map(map>0));
    x(k, 16) = max(map(map>0));
    x(k, 17) = mean(map(map>0));
    x(k, 18) = std(map(map>0));
  end
  headers= ['file, 2.5th percentile, 5th percentile, 10th percentile,25th percentile, 30th percentile, 35th  percentile, 50th percentile, 65th percentile, 70th percentile, 75th percentile, 90th percentile, 95th percentile, 97.5th percentile, min, max, mean, standard deviation'];
  fid = fopen(filename,'W+');
  fprintf(fid,'%s\n', headers);
  fclose(fid);
  dlmwrite(filename, x,'-append')
end
toc
