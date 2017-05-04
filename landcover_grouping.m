% Clear all data in previous running and close all figures.
clear; close all; clc;

%%% landcover_grouping.m : grouping the classes of land cover.
%% Input : tiff images labeled by classes, [1:16]
%% Output: tiff images grouping the classes, [1:4]
%% Classes mapping rules:
%%   forests     1:5   --> 1
%%   savannas    6:9   --> 2
%%   grasslands  10    --> 3
%%   croplands   12&14 --> 4
%%   others      [0, 11, 13, 15, 16] --> 0


% Define NODATA values
% NODATA_LANDCOVER =  -128;

years = 2001:2013
basename = '/home2/dongmeic/fire/data/LC_China/';
outputfolder = '/home2/dongmeic/fire/output/LC_China/';

% Make the output dir if it does not exist.
if ~exist(outputfolder)
  mkdir(outputfolder);
end

for y=years
  %% Process the image in every year.
  filename = [basename 'LC_China_' num2str(y) '.tif'];
  [landcover ref] = geotiffread(filename);
  output = landcover;

  %% Mapping classes into groups
  output(landcover == 1)  = 1;
  output(landcover == 2)  = 1;
  output(landcover == 3)  = 1;
  output(landcover == 4)  = 1;
  output(landcover == 5)  = 1;
  output(landcover == 6)  = 2;
  output(landcover == 7)  = 2;
  output(landcover == 8)  = 2;
  output(landcover == 9)  = 2;
  output(landcover == 10) = 3;
  output(landcover == 12) = 4;
  output(landcover == 14) = 4;
  output(landcover == 0)  = 0;
  output(landcover == 11) = 0;
  output(landcover == 13) = 0;
  output(landcover == 15) = 0;
  output(landcover == 16) = 0;
% output(landcover == NODATA_LANDCOVER) = NODATA_LANDCOVER;
  filename = [outputfolder 'LC_China_recl_' num2str(y) '.tif'];
  geotiffwrite(filename, output, ref);
end