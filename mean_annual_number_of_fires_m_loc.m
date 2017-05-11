close all; clear; clc;
tic
%%% mean_annual_number_of_fires_m.m : calculate mean annual number of fires and Gini index
%% Gini index: see http://www.fao.org/docs/up/easypol/329/gini_index_040en.pdf page 9 function 2
%% Flood-Fill algorithm reference:  Archibald & Roy 2009, Identifying individual fires from satellite-derived burned area data 
%% Input : monthly burnt area maps
%% Output: 0.5 degree mean annual number of fires, inter-annual coefficient of variance in number of fires, Gini index
%% updated on 2017.05.10

%% Problem description: the output of Gini index is with different length in the burned area

cutoff = 8;

years = 2001:2016;
months = 1:12;

basedir = '/Users/dongmeichen/BA_China/';
outputfolder = '/Users/dongmeichen/output/results/';

%% Region of Interest
[ROI ref] = geotiffread('/Users/dongmeichen/masknew.tif');
[m n] = size(ROI);


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

number_of_fires_y = zeros(mm, nn, length(years)); % number of fires in one year
number_of_months_fire_y = zeros(mm, nn, length(years)); % number of fires in one month
fireshist = cell(mm, nn);

for idx = 1:length(years)
  y=years(idx);
  number_of_fires_m = zeros(mm, nn, length(months));
  number_of_months_fire_m = zeros(mm, nn, length(months));

  for jdx = 1:length(months)
    v =months(jdx);

    map = geotiffread([basedir 'BA_China_' num2str(y) '-' num2str(v) '.tif']);
    fires= zeros(m, n);

    %% Find all the fires points.
    [I J] = find((map > 0) & ROI);
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
	    if (map(ii, jj) > 0) && ROI(ii,jj) == 1
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
    fprintf('Got the %d year %d month fires!\n', y, v);
 
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
        number_of_fires_m(ii, jj, jdx) = c;
	if c > 0
	  number_of_months_fire_m(ii, jj, jdx) = 1;
	end
        fireshist{ii, jj} = [fireshist{ii, jj} d];
      end
    end
    fprintf('Got the %d year %d month number of fires!\n', y, v);
  end
  number_of_fires_y(:,:,idx) = sum(number_of_fires_m, 3);
  number_of_months_fire_y(:,:,idx) = sum(number_of_months_fire_m, 3);
end
toc

tic
%% Gini index
Gini = -ones(mm, nn);
firesize = -ones(mm,nn);
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
      f = b*0.25;
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
      firesize(ii,jj) = prctile(f,95)*0.25;
    end
  end
end
fprintf('Gini index finished!\n');
toc

tic
manof = -ones(mm, nn); % mean annual number of fires
anofstd = -ones(mm, nn); % standard deviation of annual number of fires
anofcv  = -ones(mm, nn); % inter-annual coefficient of variation in number of fires
nofmf = -ones(mm, nn); 
manof = mean(number_of_fires_y, 3);
anofstd = std(number_of_fires_y, 0, 3);
nofmf = sum(number_of_months_fire_y, 3);
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
filename4 = [outputfolder 'number_of_months_fire.tif'];
filename5 = [outputfolder 'maximum_fire_size.tif'];

geotiffwrite(filename1, manof, outputref);
geotiffwrite(filename2, Gini, outputref);
geotiffwrite(filename3, anofcv, outputref);
geotiffwrite(filename4, nofmf, outputref);
geotiffwrite(filename5, firesize, outputref);
toc
