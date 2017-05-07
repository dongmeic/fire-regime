close all; clear; clc;

%%% percentage_landcover_affected_m.m : calculate percentage landcover affected by fire based on monthly active fire data
%% Input : monthly burnt area maps; yearly landcover reclassfied maps
%% Output: 0.5 degree percentage landcover affected by fire
%% Description: 
%% Created by: Dongmei CHEN, 2017.05.05

tic

years = 2001:2013;
months = 1:12;
baseDir_AF   = '/home2/dongmeic/fire/output/AF_China/';
baseDir_LC   = '/home2/dongmeic/fire/output/LC_China/';
outputFolder = '/home2/dongmeic/fire/output/results/';
%% Nodata = 32767;
[ROI ref] = geotiffread([baseDir_AF 'mask.tif']);
[m n] = size(ROI);


%%% Block matrix of 100-by-100 grids.
%% Indexes of grid points, two rows:
%%   row 1 are starting points
%%   row 2 are ending points.
xIndex = [[1:100:m]; [1:100:m]+99];
xIndex(end) = m;
yIndex = [[1:100:n]; [1:100:n]+99];
yIndex(end) = n;
mm = length(xIndex); %% mm, nn, ii, jj are the indexes of grid windows.
nn = length(yIndex); %% Well, m, n, i, j are the indexes of pixels.

outputref = ref;
outputref.RasterSize = [mm nn];

lcover_yearly = zeros(1, 4, m, n);

for y=years 
  for v=months
    fprintf('Processing the %d year %d month!\n', y, v);
    burntArea = geotiffread([baseDir_AF 'afd' num2str(y) '_' num2str(v) '.tif']);
    %% burntArea(burntArea == Nodata) = 0; 
    landCover = geotiffread([baseDir_LC 'LC_China_recl_' num2str(y) '.tif']);

    for jj = 1:nn
      for ii = 1:mm
        lcover_grid = zeros(1,4);
        i = xIndex(:, ii);
        j = yIndex(:, jj);
        burntArea_grid = burntArea(i(1):i(2), j(1):j(2));
        landCover_grid = landCover(i(1):i(2), j(1):j(2));
        ROI_grid = ROI(i(1):i(2), j(1):j(2));
        x = [];
        for class = 1:4
          lcover_grid(class) = numel(find(landCover_grid((burntArea_grid > 0) & ROI_grid) == class));
        end
        lcover_yearly(:, :, i(1), j(1)) = lcover_yearly(:, :, i(1), j(1)) + lcover_grid;
      end
    end
  end
end
toc

tic
%% Find most affected landcover type
pctf = -ones(mm, nn);
pctsa = -ones(mm, nn);
pctg = -ones(mm, nn);
pctc = -ones(mm, nn); 
for jj = 1:nn
  for ii = 1:mm
    i = xIndex(:, ii);
    j = yIndex(:, jj);

    lcover_grid = lcover_yearly(:, :, i(1), j(1));
    if any(lcover_grid)
      forest = (lcover_grid(1)/sum(lcover_grid))*100;
      savanna = (lcover_grid(2)/sum(lcover_grid))*100;
      grassland = (lcover_grid(3)/sum(lcover_grid))*100;
      cropland = (lcover_grid(4)/sum(lcover_grid))*100;
      pctf(ii, jj) = forest;
      pctsa(ii, jj) = savanna;
      pctg(ii, jj) = grassland;
      pctc(ii, jj) = cropland;
    end
  end
end
toc

tic
%% Output
fprintf('Start writing output!\n');
% pctf = uint16(pctf);
% pctsa = uint16(pctsa);
% pctg = uint16(pctg);
% pctc = uint16(pctc);

filename1 = [outputFolder 'percentage_forest_affected.tif'];
filename2 = [outputFolder 'percentage_savanna_affected.tif'];
filename3 = [outputFolder 'percentage_grassland_affected.tif'];
filename4 = [outputFolder 'percentage_cropland_affected.tif'];

geotiffwrite(filename1, pctf, outputref);
geotiffwrite(filename2, pctsa, outputref);
geotiffwrite(filename3, pctg, outputref);
geotiffwrite(filename4, pctc, outputref);

toc

