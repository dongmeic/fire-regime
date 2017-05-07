close all; clear; clc;
tic
%%% mean_annual_number_of_fires.m : calculate mean annual number of fires
%% Input : yearly burnt area maps
%% Output: 0.5 degree mean annual number of fires
%% Description:
%% Created by: Peng LIN
%% Modified by: Dongmei CHEN

%% NODATA for burnt area
NODATA_BA = 65535;
cutoff = 8;

years = 2001:2016;

basedir = '/media/jane/TOSHIBA EXT/Fire/Fire Data/MODIS_BurntArea/MCD64A1/yearly/mask/';
outputfolder = '/media/jane/TOSHIBA EXT/Fire/Fire Data/MODIS_BurntArea/MCD64A1/output/';

%% basedir = './yearly_burnt_area/mask/';
%% outputfolder = './yearly_burnt_area/output/';

[map ref] = geotiffread([basedir '2001.tif']);
[m n] = size(map);
ROI = (map ~= NODATA_BA);

neighbor = [-1 -1; 0 -1; +1 -1; -1 0; +1 0; -1 +1; 0 +1; +1 +1];

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

number_of_fires_all = zeros(mm, nn, length(years));
number_of_years_fire = zeros(mm, nn, length(years));
fireshist = cell(mm, nn);

for idx = 1:length(years)
  y=years(idx);
  map = geotiffread([basedir num2str(y) '.tif']);
  fires= zeros(m, n);

  %% Find all the fires points.
  [I J] = find((map ~= 0) & ROI);
  len = length(I);
  %% Flood-Fill algorithm
  fireID = 0;
  for k = 1:len
    i = I(k);
    j = J(k);
    if fires(i, j) == 0
      fireID = fireID + 1;
      fires(i, j) = fireID;
      queue = [i, j];
      while queue
        i = queue(1, 1);
        j = queue(1, 2);
        for kk = 1:8
	  ii = i + neighbor(kk, 1);
	  jj = j + neighbor(kk, 2);
	  if (map(ii, jj) ~= 0) && (map(ii, jj) ~= NODATA_BA)
	    %% Be careful the type of map is uint16. Should turn it to double. See class(map);
	    if (abs(double(map(i, j)) - double(map(ii, jj))) <= cutoff) && (fires(ii, jj) == 0)
	      fires(ii, jj) = fires(i, j);
	      queue = [queue; ii jj];
	    end
	  end
        end
        queue(1, :) = [];
      end
    end
  end
  fprintf('Got the %d yearly fires!\n', y);

  number_of_fires = zeros(mm, nn);
  fires_grid = [];
  d = [];

  for jj = 1:nn
    for ii = 1:mm
      i = xIndex(:, ii);
      j = yIndex(:, jj);
      fires_grid = fires(i(1):i(2), j(1):j(2));
      a = fires_grid(fires_grid > 0);
      b = unique(a);
      c = numel(b);
      d = zeros(1, c);
      for e = 1:c
	    d(e) = sum(a == b(e));
      end
      number_of_fires_all(ii, jj, idx) = c;
	  if c > 0
	    number_of_years_fire(ii, jj, idx) = 1;
	  end
      fireshist{ii, jj} = [fireshist{ii, jj} d];
    end
  end
  fprintf('Got the %d year number of fires!\n', y);
end
toc

tic
%% Gini index
Gini = zeros(mm, nn);
a = [];

for jj = 1:nn
  for ii = 1:mm
    i = xIndex(:, ii);
    j = yIndex(:, jj);
    if any(fireshist{ii, jj})
      a = fireshist{ii, jj};
      b = sort(a);
      c = sum(b);
      e = numel(b);
      d = zeros(1,e);
      q = zeros(1,e);
      p = zeros(1,e);
      z = zeros(1,e);
      d(1)=b(1)/c;
      q(1)=d(1);
      p(1)=1/e;
      z(1)=q(1)*p(1);
      for k=2:e
        d(k)=b(k)/c;
        q(k)=d(k)+q(k-1);
        p(k)=k/e;
        z(k)=(q(k)+q(k-1))*(p(k)-p(k-1));
      end
      Gini(ii, jj) = 1-sum(z);
    end
  end
end
fprintf('Gini index finished!\n');
toc

tic
manof = zeros(mm, nn);
anofstd = zeros(mm, nn);
anofcv  = zeros(mm, nn);
nofyf = zeros(mm, nn);
manof = mean(number_of_fires_all, 3);
anofstd = std(number_of_fires_all, 0, 3);
nofyf = sum(number_of_years_fire, 3);
nof = sum(number_of_fires_all, 3); 

%% Coefficient of variance
for jj = 1:nn
  for ii = 1:mm
    if manof(ii, jj) ~= 0
      anofcv(ii, jj) = anofstd(ii, jj) / manof(ii, jj);
    end
  end
end
fprintf('Got coefficient of variance!\n');
toc

tic
%% Output
filename1 = [outputfolder 'mean_annual_number_of_fires.tif'];
filename2 = [outputfolder 'Gini_index.tif'];
filename3 = [outputfolder 'mean_annual_number_of_fires_cv.tif'];
filename4 = [outputfolder 'number_of_years_fire.tif'];
filename5 = [outputfolder 'number_of_fires.tif'];

% geotiffwrite(filename1, manof, outputref);
% geotiffwrite(filename2, Gini, outputref);
% geotiffwrite(filename3, anofcv, outputref);
% geotiffwrite(filename4, nofyf, outputref);
geotiffwrite(filename5, nof, outputref);
toc
