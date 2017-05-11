
  y=years(idx);
  number_of_fires_m = -ones(mm, nn, length(months));
  number_of_months_fire_m = -ones(mm, nn, length(months));

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
	    if (ii>0 && jj>0 && ii<=m && jj <=n) && (map(ii, jj) > 0) && ROI(ii,jj) == 1
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
    fires1 = fires;
    
    
    
    %% Flood-Fill algorithm
    fires= zeros(m, n);
    fireID = 0;
    for k = 1:len
      i = I(k);
      j = J(k);
      if fires(i, j) == 0
        fireID = fireID + 1;
        fires(i, j) = fireID;
        queue = [i, j];
        
          i = queue(1, 1);
          j = queue(1, 2);
          for kk = 1:8
	    ii = i + neighbor(kk, 1);
	    jj = j + neighbor(kk, 2);
	    if (ii>0 && jj>0 && ii<=m && jj <=n) && (map(ii, jj) > 0) && ROI(ii,jj) == 1
	      %% Be careful the type of map is uint16. Should turn it to double. See class(map);
	      if (abs(double(map(i, j)) - double(map(ii, jj))) <= cutoff) && (fires(ii, jj) == 0)
	        fires(ii, jj) = fires(i, j);
	        queue = [queue; ii jj];
	      end
	    end
          end
        end
    end
    fires2 = fires;

    
    
    
    %% Flood-Fill algorithm
    fires= zeros(m, n);
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
	    if (ii>0 && jj>0 && ii<=m && jj <=n) && (map(ii, jj) > 0) && ROI(ii,jj) == 1
	      %% Be careful the type of map is uint16. Should turn it to double. See class(map);
	      if (abs(double(map(i, j)) - double(map(ii, jj))) <= cutoff) && (fires(ii, jj) == 0)
	        fires(ii, jj) = fires(i, j);
	        queue = [queue; ii jj];
	      end
	    end
          end
          queue = queue(2:end,:);
        end
      end
    end
    fires3 = fires;

sum( fires1(:) > fires2(:) | fires1(:) < fires2(:) )
sum( fires1(:) > fires3(:) | fires1(:) < fires3(:) )
