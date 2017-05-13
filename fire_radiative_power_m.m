close all; clear; clc;
tic
%%% fire_radiative_power_m.m : calculate 95th frp
%% Input : yearly active fire maps
%% Output: 0.5 degree 95th percentile of fire radiative power
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
    basedir2 = '/home2/dongmeic/fire/';
    outputfolder = '/home2/dongmeic/fire/output/results/';
  end
end

years = 2001:2015;
months=1:12;

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

frp_all = cell(mm, nn);

for idx = 1:length(years)
	y=years(idx);
	for jdx = 1:length(months)
		v=months(jdx);
		fprintf('Processing activefire the %d year %d month!\n', y,v);
		map = geotiffread([basedir 'afd' num2str(y) '_' num2str(v) '_frp.tif']);
		frp = cell(mm, nn);
		for jj = 1:nn
			for ii = 1:mm
        if (ii>0 && jj>0 && ii<=m && jj <=n) && ROI(ii,jj) == 1
          i = xIndex(:, ii);
          j = yIndex(:, jj);
          map_grid = map(i(1):i(2), j(1):j(2));
          ROI_grid = ROI(i(1):i(2), j(1):j(2));
          map_grid(ROI_grid == 0) = 0;
          frp{ii, jj} = map_grid(map_grid > 0);
          frp_all{ii,jj} = [frp_all{ii,jj};frp{ii, jj}];
        end
			end
		end
	end
  fprintf('Got the %d year frp!\n', y);
end
toc

tic
frp_pct = -ones(mm, nn);
for jj = 1:nn
  for ii = 1:mm
    if any(frp_all{ii,jj})  % If this grid has any fire.
      frp_pct(ii,jj) = prctile(frp_all{ii,jj},95);
    end
  end
end
toc

tic
%% Output
fprintf('Start writing output!\n');
filename = [outputfolder 'fire_radiative_power.tif'];
geotiffwrite(filename, frp_pct, outputref);
toc
