tic
close all; clear; clc;

%%% mean_annual_active_fire_m.m : calculate mean annual area burned
%% Input : monthly burnt area maps
%% Output: 0.5 degree mean annual burn area map
%% Description:
%% Created by: Dongmei CHEN, 2017.05.05
%% updated on 2017.05.10

basedir = '/Users/dongmeichen/BA_China/';
outputfolder = '/Users/dongmeichen/output/results/';
%% Nodata = 32767;

years = 2001:2016;
months = 1:12;

%% Map size and reference system
[ROI ref] = geotiffread('/Users/dongmeichen/masknew.tif');
[m n] = size(ROI);

xIndex = [[1:100:m]; [1:100:m]+99];
xIndex(end) = m;
yIndex = [[1:100:n]; [1:100:n]+99];
yIndex(end) = n;
mm = length(xIndex); %% mm, nn, ii, jj are the indexes of grid windows.
nn = length(yIndex); %% Well, m, n, i, j are the indexes of pixels.

outputref = ref;
outputref.RasterSize = [mm nn];

afstd = zeros(mm, nn);
afcv = -ones(mm, nn);
afmean = zeros(mm, nn);
gridAll = zeros(mm, nn);
afAll = zeros(mm, nn);
afpct = zeros(mm, nn);
burnPixel_y = zeros(mm, nn, length(years));
gridSize = zeros(mm, nn, length(years));
% gridSize_y = zeros(mm, nn, length(years));
%% Extract the data
for idx = 1:length(years)
  y=years(idx);
  burnPixel_m = zeros(mm, nn, length(months));
  % gridSize_m = zeros(mm, nn, length(months));
  for jdx = 1:length(months)
    v=months(jdx);
    fprintf('Processing burned area the %d year %d month!\n', y, v);
    map = geotiffread([basedir 'BA_China_' num2str(y) '-' num2str(v) '.tif']);
    %% map(map == Nodata) = 0; 
    for jj = 1:nn
      for ii = 1:mm
        i = xIndex(:, ii);
        j = yIndex(:, jj);
        map_grid = map(i(1):i(2), j(1):j(2));
        ROI_grid = ROI(i(1):i(2), j(1):j(2));
        map_grid(ROI_grid == 0) = 0;
        if sum(sum(map_grid > 0))>0
            burnPixel_m(ii, jj, jdx) = sum(sum(map_grid > 0));
        end
        % gridSize_m(ii, jj, jdx) = (i(2)-i(1)+1)*(j(2)-j(1)+1);
        gridSize(ii, jj, idx) = (i(2)-i(1)+1)*(j(2)-j(1)+1);
      end
    end
  end
  burnPixel_y(:,:,idx) = sum(burnPixel_m, 3);
  % gridsize_y(:,:,idx) = sum(gridSize_m, 3);
  fprintf('Burned area the %d year finished!\n', y);
end


afmean = mean(burnPixel_y, 3);
afAll = sum(burnPixel_y, 3);
gridAll = sum(gridSize, 3);
afstd = std(burnPixel_y, 0, 3);

for jj = 1:nn
  for ii = 1:mm
    afpct(ii, jj) = (afAll(ii, jj) / gridAll(ii, jj))*100; 
    if afmean(ii, jj) > 0
      afcv(ii, jj) = afstd(ii, jj) / afmean(ii, jj);
    end
  end
end

%afpct = uint16(afpct);
filename1 = [outputfolder 'mean_annual_area_burned.tif'];
geotiffwrite(filename1, afmean, outputref);
filename2 = [outputfolder 'mean_annual_area_burned_cv.tif'];
geotiffwrite(filename2, afcv, outputref);
filename3 = [outputfolder 'pixels_burned_percentage_BA.tif'];
geotiffwrite(filename3, afpct, outputref);

toc
