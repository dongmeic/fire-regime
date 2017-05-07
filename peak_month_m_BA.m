close all; clear; clc;
tic

% peak_month_m: calculate the month when the maxium burned pixels occurred by adding all burned pixels month by month
%% Output: 0.5 degree peak month and maxium burned pixels in one month
%% Created by: Dongmei CHEN, 2017.05.05

%% Region of Interest
basedir = '/home2/dongmeic/fire/output/BA_China/';
outputfolder = '/home2/dongmeic/fire/output/results/';

years = 2001:2016;
months = 1:12;

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

pixels_years = zeros(mm, nn, length(years));
pixels_months = zeros(mm, nn, length(months));
for idx = 1:length(months)
  v = months(idx);
  for jdx = 1:length(years)
    y = years(jdx);
    fprintf('Processing the %d year %d month!\n', y, v);
    map = geotiffread([basedir 'BA_China_' num2str(y) '-' num2str(v) '.tif']);
    %% map(map == Nodata) = 0; 
    for jj = 1:nn
      for ii = 1:mm
        i = xIndex(:, ii);
        j = yIndex(:, jj);
        map_grid = map(i(1):i(2), j(1):j(2));
        ROI_grid = ROI(i(1):i(2), j(1):j(2));
        map_grid(ROI_grid == 0) = 0;
        pixels_years(ii, jj, jdx) = sum(sum(map_grid > 0));
      end
    end
  end
  pixels_months(:, :, idx) = sum(pixels_years, 3);
  fprintf('Got the %d month all years area burned!\n', v);
end
toc

tic
peak = zeros(mm, nn);
month = zeros(mm, nn);
[peak, month] = max(pixels_months, [], 3);

% Check the number of peaks
if 0
  peaks_number = zeros(mm, nn);
  for k=1:months
    peaks_number = peaks_number + (pixels_months(:, :, k) == peak);
  end
  fprintf(peaks_number);
end

% Length of fire period(lfp)
firepixels = zeros(mm, nn);
firepixels = sum(pixels_months,3);
sf = zeros(mm, nn);
lfp = -ones(mm,nn);

for jj = 1:nn
    for ii = 1:mm       
        if firepixels(ii, jj) ~= 0 
            for i=1:12
                if pixels_months(ii, jj, i) >= 0.1 * firepixels(ii, jj)
                    sf(ii,jj)=sf(ii,jj)+1;
                end                  
            end
            lfp(ii,jj)= sf(ii,jj);
        end
    end
end

% Output
filename1 = [outputfolder 'peak_month_burned_pixels_BA.tif'];
geotiffwrite(filename1, peak, outputref);
peaks_ROI = (peak ~= 0);
month(~peaks_ROI)= 0;
filename2 = [outputfolder 'peak_month_BA.tif'];
geotiffwrite(filename2, month, outputref);
filename3 = [outputfolder 'length_of_fire_period_BA.tif'];
geotiffwrite(filename3, lfp, outputref);
toc


