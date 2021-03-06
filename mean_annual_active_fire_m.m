tic
close all; clear; clc;

%%% mean_annual_active_fire_m.m : calculate mean annual area burned
%% Input : monthly burnt area maps
%% Output: 0.5 degree mean annual burn area map
%% Description:
%% Created by: Dongmei CHEN, 2017.05.05

% working in different computers without having to manually change the paths
if isequal( getenv('UserName') , 'Andrea' )
  basedir = '';
  basedir2 = '';
  outputfolder = 'results/';
else
  [~,computer_name]=system('hostname');
  if length(computer_name)>=20 && isequal( computer_name(1:20), 'd136-228.uoregon.edu' )
    basedir = '/Volumes/dongmeichen/output/AF_China/';
    basedir2 = '/Volumes/dongmeichen/';
    outputfolder = '/Volumes/dongmeichen/output/results/';
  else
    basedir = '/home2/dongmeic/fire/output/AF_China/';
    outputfolder = '/home2/dongmeic/fire/output/results/';
    basedir2 = '/home2/dongmeic/fire/';
  end
end


%% Nodata = 32767;

years = 2001:2015;
months = 1:12;

%% Map size and reference system
[ROI ref] = geotiffread([basedir2 'masknew.tif']);
ROI( ROI < 0 ) = 0;
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
    fprintf('Processing active fire the %d year %d month!\n', y, v);
    map = geotiffread([basedir 'afd' num2str(y) '_' num2str(v) '.tif']);
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
  fprintf('Active fire the %d year finished!\n', y);
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
filename1 = [outputfolder 'mean_annual_active_fire.tif'];
geotiffwrite(filename1, afmean, outputref);
filename2 = [outputfolder 'mean_annual_active_fire_cv.tif'];
geotiffwrite(filename2, afcv, outputref);
filename3 = [outputfolder 'pixels_burned_percentage.tif'];
geotiffwrite(filename3, afpct, outputref);

toc
