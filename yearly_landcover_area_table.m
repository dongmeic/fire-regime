% Yearl land cover area extracted from MODIS product MCD12Q1

% Define NODATA values
NODATA_LANDCOVER =  -128;

years = 2001:2013
lcover_statistics = [];

% Extract Data
for y=years
  
  landcover = geotiffread(strcat('/Volumes/dongmeic/fire/output/LC_China/LC_China_recl_', num2str(y), '.tif'));  

  % Region of interest
  ROI = (landcover ~= NODATA_LANDCOVER); % when map is not equal to NODATA value.
  
  % land cover statistics
  x = [];
  for class = 0:4
    x(class+1,1) = class;
    x(class+1,2) = numel(find(landcover(ROI) == class)); 
  end
  
  lcover_statistics = [lcover_statistics x(:,2)*25]; % feeding the output, number of pixels burnt for each land cover class
  fprintf('%d processed!\n', y);
end
lc_table=[x(:,1)';lcover_statistics'];
lc_table=[[0;years'] lc_table];
csvwrite('/Volumes/dongmeic/fire/output/revision/landcover_area_table.csv',lc_table)
