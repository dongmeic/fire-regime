close all; clear; clc;
tic
%%% fire_season_duration.m : calculate fire season duration, which is defined as the period 
%% between the day when 10% of fire pixels in total occurred and the day when 90% of fire pixels accumulated; 
%% peak is the week when maximum fire pixels accumulated in one week
%% Input : monthly burnt area maps
%% Output: 0.5 degree fire season duration

%% Problem description: the output of fire season duration is with different length in the burned area 
%% from the previous output (e.g., burned area, number of fires) and the boxplot of data is not clear


basedir = '/home2/dongmeic/fire/output/BA_China/';
outputfolder = '/home2/dongmeic/fire/output/results/';

years = 2001:2016;
months = 1:12;

%% Region of Interest
[ROI ref] = geotiffread([basedir 'mask.tif']);
[m n] = size(ROI);

xIndex = [[1:100:m]; [1:100:m]+99];
xIndex(end) = m;
yIndex = [[1:100:n]; [1:100:n]+99];
yIndex(end) = n;
mm = length(xIndex); %% mm, nn, ii, jj are the indexes of grid windows.
nn = length(yIndex); %% Well, m, n, i, j are the indexes of pixels.

outputref = ref;
outputref.RasterSize = [mm nn];

dayshist_all = zeros(1, 366, mm, nn);

%% reading data (importing them from files)
for y=years
  for v=months
    fprintf('Processing the %d year %d month!\n', y, v);
    map = geotiffread([basedir 'BA_China_' num2str(y) '-' num2str(v) '.tif']);
    for jj = 1:nn
      for ii = 1:mm
        i = xIndex(:, ii);
        j = yIndex(:, jj);
        dayshist= zeros(1, 366);
        map_grid = map(i(1):i(2), j(1):j(2));
        ROI_grid = ROI(i(1):i(2), j(1):j(2));
        [I J] = find((map_grid > 0) & ROI_grid);
        len = length(I);
        for k = 1:len
          a = I(k);
          b = J(k);
          day = map_grid(a, b);
          dayshist(day) = dayshist(day) + 1;
        end
        dayshist_all(:, :, ii, jj) = dayshist_all(:, :, ii, jj) + dayshist; 
      end
    end
  end
  fprintf('Got the %d year dayshist!\n', y);
end
toc

tic
burntdates = zeros(1, 366);
fs_duration = zeros(mm, nn);
pw_start = zeros(mm, nn);
pw_number = zeros(mm, nn);
number_of_fs_duration_equal_to_0 = 0;
for jj = 1:nn
  for ii = 1:mm
    burntdates = dayshist_all(:, :, ii, jj);
    if any(burntdates)  % If this grid has any fire.
      for l = 2:366
        burntdates(l) = burntdates(l) + burntdates(l - 1);
      end

      a10 = burntdates(366) * 0.10;
      a90 = burntdates(366) * 0.90;
      day10 = 1;
      while burntdates(day10) < a10
        day10 = day10 + 1;
      end
      day90 = 366;
      while burntdates(day90) > a90
        day90 = day90 - 1;
      end
      fs_duration(ii, jj) = day90 - day10;
      if fs_duration(ii, jj) == 0
        number_of_fs_duration_equal_to_0 = number_of_fs_duration_equal_to_0 + 1;
      end
      if fs_duration(ii, jj) == 365
        fprintf('Error: no fires at all: ii = %d, jj = %d\n', ii, jj)
        return;
      end
      

      p = day10;
      dayw = [];
      while p + 7 <= day90
        w = burntdates(p+7)-burntdates(p);
        p = p + 7;
        dayw = [dayw w];
      end
      if any(dayw)
        [peak peak_start] = max(dayw);
        pw_start(ii, jj) = peak_start;
        pw_number(ii,jj) = sum(dayw == peak);
      end
    end
  end
end
toc

%% Andrea: check non-zero and non-empty points in several variables...
dayshist_all_greater_than_0 = zeros(mm, nn);
for jj = 1:nn
	for ii = 1:mm
		if any( dayshist_all(:, :, ii, jj) )
			dayshist_all_greater_than_0( ii, jj ) = 1;
		end
	end
end
disp( ['Fire season duration is 0 in ' num2str(number_of_fs_duration_equal_to_0) ' places where you do have fire measurements!'] )
disp( ['Total number of places where fire season duration is different from 0: ' num2str( sum(fs_duration(:)~=0) ) ] )
disp( ['Total number of places where dayshist is non-empty: ' num2str( sum( dayshist_all_greater_than_0(:) ))])

tic
%% Output
fprintf('Start writing output!\n');
fs_duration = uint16(fs_duration);
pw_start = uint16(pw_start);
pw_number = uint16(pw_number);

filename1 = [outputfolder 'fire_season_duration.tif'];
filename2 = [outputfolder 'peak_start.tif'];
filename3 = [outputfolder 'number_of_peaks.tif'];
geotiffwrite(filename1, fs_duration, outputref);
geotiffwrite(filename2, pw_start, outputref);
geotiffwrite(filename3, pw_number, outputref);
toc
