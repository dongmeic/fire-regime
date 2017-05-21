close all; clear; clc;

tic

savemyresults = false;

%%% kmeans_clustering.m : clustering by applying kmeans
%% Input: spatial,temporal and intensity attributes of fire regimes
%% Output: cluster results
%% Description: need function imenlarger
%% Created by: Peng LIN
%% Modified by: Dongmei CHEN, Andrea MASIERO

basedir = '/Volumes/dongmeic-10/fire/output/revision/results/';

%% Pre-defined values.
nCluster = 6; % decided by 'elbow' method
nPCA = 12; % without PCA
%% File names.
filenames = {
	     'mean_annual_number_of_fires.tif'; 
         	    %'mean_annual_area_burned.tif'; 
	     %'mean_annual_area_burned_cv.tif'; 
                      %'maximum_fire_size.tif';
          'mean_annual_active_fire.tif';
          'mean_annual_number_of_fires_cv.tif'; 
	     'mean_annual_active_fire_cv.tif'; 
	     'fire_season_duration.tif'; 
         %'length_of_fire_period.tif';
         %'length_of_fire_period_BA.tif';
             %'peak_start.tif'; 
             'peak_month.tif';
             'fire_radiative_power.tif'; 
              'Gini_index.tif';
	     'percentage_forest_affected.tif';
	     'percentage_savanna_affected.tif';
	     'percentage_grassland_affected.tif';
	     'percentage_cropland_affected.tif';
        %'percentage_forest_affected_BA.tif';
       %'percentage_savanna_affected_BA.tif';
       %'percentage_grassland_affected_BA.tif';
       %'percentage_cropland_affected_BA.tif';       
}

file_name = {
	     'mean annual number of fires';
         %'mean annual burned area';
	     %'CV of annual burned area';
                %'maximum fire size';
          'mean annual active fire';
         'CV of annual number of fires';
	     'CV of annual active fire';
	     'fire season duration';
         %'length of fire period';
		     %'peak start';
         'peak month';
         'fire radiative power';
         	     'Gini index';
	     'percentage forest affected';
	     'percentage savanna affected';
	     'percentage grassland affected';
	     'percentage cropland affected';	     
}

[map ref] = geotiffread([basedir filenames{1}]);
[m n] = size(map);
len = length(filenames);

InputData = zeros(m * n, len);
size(InputData)

for k = 1:len
  map = geotiffread([basedir filenames{k}]);
  if ~all(size(map) == [m n])
    fprintf('[Error]: size is mismatched!\n');
  end
  InputData(:, k) = map(:);
%   mincut=prctile(InputData(:, k),2);
%   maxcut=prctile(InputData(:, k),98);
%   tmp = InputData(:, k);
%   tmp(tmp <= mincut) = 0;
%   tmp(tmp >= maxcut) = 0;
%   InputData(:, k) = tmp;

% set no data value for variables(mean annual number of fires (area burned), mean annual
% active fire, fire season duration and peak month)
  if k==1 || k==2 || k==5 || k==6
      tmp = InputData(:, k);
      tmp(tmp==0)=-1;
      InputData(:, k) = tmp;
  end
end


%% reset the input
X = InputData;
% exclude samples with no data value in all variables 
ind = all( X(:,:)==-1, 2 );
X = X( ~ind , : );
ind = ~ind;


%% rescale the data
% logarithm
% X = log(1+X);
type_of_normalization = 3;
switch type_of_normalization
 case 1
% divided by the maximum value suggested by Andrea
  for k = 1:len
    X(:, k) = X(:, k)/max(max(X(:, k)));
  end
  % mu = mean(X);
  % Xn = X - repmat(mu, size(X,1), 1);
  Xn = X;

% by Andrea
 case 2
  mu = mean(X);
  Xn = X - repmat(mu, size(X,1), 1);
  Xn = Xn ./ repmat( std(Xn) , size(X,1), 1  );

 case 3
  % another way to normalize data
  for k = 1:len
    X(:, k) = (X(:, k)-min(min(X(:, k))))/(max(max(X(:, k)))-min(min(X(:, k))));
  end
  %mu = mean(X);
  %Xn = X - repmat(mu, size(X,1), 1);
  Xn = X;
end
  

if 0
	%look at the histogram of normalized data
	for i = 1:len
		figure(132)
		clf
		hist(X(:,i))
		if i < len
			title(['Histgram: ' file_name{i}]) % press any button
			pause
		end
	end
end

% See 'eblow' to decide the number of clusters
if savemyresults
    kcluster = 20;
    wss = zeros(1, kcluster);
    for k = 1:kcluster
        [idx,C,sumd,D] = kmeans(Xn, k,'MaxIter',1000,'Display','final','Replicates',100);
        wss(k) = norm(sumd)^2;
        %wss(k) = sum(sumd);
    end
    plot(wss,'LineWidth', 1.5);
    s=sum(wss);
    wssn=wss/s;
    tmp = wss / sum(wss)
    ttmp = tmp(1:(end-1)) - tmp(2:end);
    fprintf('The relative differences between within-cluster distances sum of squares are following:\n');
    ttmp
    set(gca,'FontSize', 12, 'FontWeight', 'Bold');
    xlabel('Number of Clusters','FontSize', 12, 'FontWeight', 'Bold');
    ylabel('Within-cluster Distances Sum of Squares','FontSize', 12, 'FontWeight', 'Bold');
    axis([1 kcluster 0 inf])
    if savemyresults
      saveas(gcf,[basedir 'clusters/kmeans_eblow_' date '.png'],'png');
    end
end

%% PCA.
Y = 1/n * (Xn' * Xn);
[U S V] = svd(Y);
s = diag(S);

%% Display the "elbow" of the sigmas of PCA.

if 0 
  figure;
  plot(s/sum(s));
  xlabel('k-th principle component');
  ylabel('weight');
  axis([0, 15, 0, 0.01]);
  set(gca,'XTick',0:1:15);
  set(gca,'YTick',0:0.005:0.05);
  saveas(gcf,[basedir 'clusters/elbow.png'],'png')
end

X_pca = Xn * V;

k_min = len;
cumu = zeros(len, 1);
for k = 1:len
  cumu(k) = sum(s(1:k)) / sum(s);
end

% Andrea: the number of principal components to be considered can be chosen
% following different types of considerations.. for instance the percentual
% error as you were doing here, or looking at how s values decreases.
% here s
if 0
  figure; 
  plot(s, 'LineWidth', 2);
  xlabel('Principle component index', 'FontSize', 12, 'FontWeight', 'Bold');
  ylabel('Weight', 'FontSize', 12, 'FontWeight', 'Bold');
  set(gca, 'FontSize', 12, 'FontWeight', 'Bold' );
  if savemyresults
	  saveas(gcf,[basedir 'clusters/elbow_' date '.png'],'png')
  end
end

% display and save PCA results
table = zeros(len, 4);
table(:, 1) = [1:len]';
table(:, 2) = s;
table(:, 3) = s / sum(s);
table(:, 4) = cumu;
fprintf('The PCA results are following:\n');
for k = 1:len
  fprintf('%2d %12.2f %1.7f %1.7f\n', table(k, 1), table(k, 2), table(k, 3), table(k, 4));
end
if savemyresults
  csvwrite([basedir 'clusters/pca_' date '.csv'], table);
end

% set C as the matrix containing the principal directions
% set XX as the matrix containing the scores
C = V(:,1:nPCA);
XX = X_pca(:,1:nPCA)';

%% computing variances of the original layer variables explained by nPCA principal components
Err = Xn - (C*XX)';
figure
plot( 100-sum( (Err).^2 ) ./ sum( Xn.^2 ) * 100, 'LineWidth', 2 )
hold on;
if nPCA <=5
  plot([1, 12], [80, 80], '--r', 'LineWidth', 1)
end
%figure,plot( sum( (C*XX).^2 , 2 )' ./ sum( Xn.^2 ) * 100 )
xlabel('Fire regime variables', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Percent of variability explained with 5 PCs', 'FontSize', 12, 'FontWeight', 'Bold');
set(gca, 'FontSize', 12, 'FontWeight', 'Bold' );
xlim([1 12])
% axis([1 12 60 100]);
if 0
  saveas(gcf,[basedir 'clusters/percent_of_variability_' date '.png'],'png');
end
%%

disp('abs(C): a high value in column i and row j means that there is a strong relation between variable j and PC i')
disp(abs(C))
%% Does it mean something with positive and negative numbers in C correlation?
if savemyresults
	csvwrite([basedir 'clusters/correlation_' date '.csv'], C);
end


%[idx, regime_classes_centers] = kmeans(Xn, nCluster);
%[idx,C,sumd,D] = kmeans(Xn, nCluster,'MaxIter',1000,'Display','final','Replicates',100);
[idx,regime_classes_centers,sumd,D] = kmeans(XX', nCluster,'MaxIter',1000,'Display','final','Replicates',100);
% colors = hot(30);
% mycolors = colors(3:3:end,:);
mycolors = [0.894 0.102 0.11; % red
  0.216 0.494 0.722; % blue
  0.302 0.686 0.29; % green
  0.596 0.306 0.639; % purple
  1 0.498 0; % orange
  1 1 0.2] % yellow 
% print a table to see the variable value range in all clusters
k=nCluster;
if savemyresults
    table = zeros(len,2*k);
    for i=1:k        
        table(:,(3*i-2):(3*i))=[min(X(idx==i,:))',median(X(idx==i,:))',max(X(idx==i,:))'];
    end
    fprintf('The variable values (min, median, max) for n clusters are following:\n');
    table
    if savemyresults
      csvwrite([basedir 'clusters/variable_values_' date '.csv'], table);
    end
end


figure(138)
clf
[silh,h] = silhouette(Xn,idx);
h = gca;
h.Children.EdgeColor = [.8 .8 1];
xlabel 'Silhouette Value';
ylabel 'Cluster';
pause
if savemyresults
    saveas(gcf,[basedir 'clusters/sihouettev' date '.png'],'png');
end

labels = {'MANF'; 'MAFD'; 'CVNF'; 'CVFD'; 'FSD'; 'FPM'; 'FRP'; 'GI'; 'PFA'; 'PSA'; 'PGA'; 'PCA'};
if savemyresults
    figure('Units','normalized')
    for i=1:len  % len = 12
        subplot(3,4,i,'replace'); %,'labels',i
        boxplot(X(:,i), idx, 'orientation','horizontal','colors',mycolors(1:k,:),'boxstyle', 'filled', 'widths',1.5, 'OutlierSize',2);
        h = findobj(gca,'tag','Median');
        set(h,'linestyle','-.');
        set(h,'Color',[0 0 0])
        xlabel([num2str(i) '-' labels{i}],'FontSize', 10, 'FontWeight', 'Bold');
%         h = findobj(gca,'Tag','Box');
%         for j=1:length(h)
%             patch(get(h(j),'XData'),get(h(j),'YData'),'y','FaceAlpha',.5);
%         end
        % boxplot(X(idx==1,6)', 'orientation','horizontal','colors',mycolors(1,:));
    end
    
    if savemyresults
      saveas(gcf,[basedir 'clusters/boxplot_cluster_variable' date '.png'],'png');
    end    
    
end

% if 0
%     figure('Units','normalized')
%     for i=1:k  % k = 7
%         subplot(4,2,i);
%         boxplot(X(idx==i,:),'colors',mycolors(i,:));
%     end
% 
%     if savemyresults
%       saveas(gcf,[basedir 'clusters/boxplot_variable_cluster' date '.png'],'png');
%     end      
% end

%see how the classes partitioned by kmeans are spatially distributed
myind = find( ind );

use_color_image = 0;
switch use_color_image
	case 0
		% grayscale image
		I1 = zeros(m,n);
		for i = 1:nCluster
			I1(myind(idx==i)) = i/nCluster;
		end
		hfig0 = figure(135);
        s0 = 6;
        imshow(imenlarger(I1,s0))
%         imshow(I1)
		if savemyresults
			print(hfig0, '-dpng', '-r500', [basedir 'clusters/regime' date '.png'])
		end
	case 1
		% color image
		
		% if you change the number of clusters you have to manually set
		% the colors that you want here
		% The i-th row of mycolors corresponds to the RGB color of the i-th
		% segmented region. RGB values are normalized to 1 (hence values
		% have to be set between 0 and 1) 
%		mycolors = [1 0 1; % magenta 101
%			1 0 0; % red 100
%			0 0 1; % blue 001
%			0 1 0; % green 010
%            0 1 1; % cyan 011
%                        1 1 0];% yellow 110
		
		I1(:,:,1) = ones(m,n); %if you want black background use zeros intead of ones
		I1(:,:,2) = ones(m,n); %if you want black background use zeros intead of ones
		I1(:,:,3) = ones(m,n); %if you want black background use zeros intead of ones
		for i = 1:nCluster
			Itmp = I1(:,:,1);
			Itmp(myind(idx==i)) = mycolors(i,1);
			I1(:,:,1) = Itmp;
			Itmp = I1(:,:,2);
			Itmp(myind(idx==i)) = mycolors(i,2);
			I1(:,:,2) = Itmp;
			Itmp = I1(:,:,3);
			Itmp(myind(idx==i)) = mycolors(i,3);
			I1(:,:,3) = Itmp;
		end
		hfig1 = figure(135);
        s1 = 6;
        imshow(imenlarger(I1, s1))
%         imshow(I1)
        %title('Fire regimes')
		if savemyresults
			print(hfig1, '-dpng', '-r500', [basedir 'clusters/regime_' date '.png'])
		end
end


if savemyresults
	% run this instruction if you want to save figure 135
	% as tiff image
	% saveas(gcf,[basedir 'clusters/regime.png'],'png');
	% print([basedir 'clusters/regimes.png'],'-dpng','-r500')
	% Andrea: unfortunately I don't have the function "geotiffwrite" in my
	% Matlab version. Anyway, if you want to save geographical information
	% in the tiff file, you need such information as well... whereas the
	% ref matrix is currently an empty matrix
	geotiffwrite([basedir 'clusters/regime_' date '.tif'], I1, ref);
end



for i = 1:nCluster
	figure(134)
	clf
	subplot(2,1,1)
	imshow(I1)
	
	switch use_color_image
		case 0
			I = zeros(m,n);
			I(myind(idx==i)) = 1;
		case 1
			I(:,:,1) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,2) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,3) = ones(m,n);%if you want black background use zeros intead of ones
			
			Itmp = I(:,:,1);
			Itmp(myind(idx==i)) = mycolors(i,1);
			I(:,:,1) = Itmp;
			Itmp = I(:,:,2);
			Itmp(myind(idx==i)) = mycolors(i,2);
			I(:,:,2) = Itmp;
			Itmp = I(:,:,3);
			Itmp(myind(idx==i)) = mycolors(i,3);
			I(:,:,3) = Itmp;
	end
	subplot(2,1,2) %% I want plot I only and save each figure in enlarger but I also want to keep this subplot when I need 
	%ANDREA: see the code some lines below (figure(143))
	imshow(I)
	if i <= nCluster
		title(['Fire regime ' num2str(i)]) % press any button
		pause
	end
end


for i = 1:nCluster
	hfig2 = figure(143);
	clf
	switch use_color_image
		case 0
			I = zeros(m,n);
			I(myind(idx==i)) = 1;
		case 1
			I(:,:,1) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,2) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,3) = ones(m,n);%if you want black background use zeros intead of ones
			
			Itmp = I(:,:,1);
			Itmp(myind(idx==i)) = mycolors(i,1);
			I(:,:,1) = Itmp;
			Itmp = I(:,:,2);
			Itmp(myind(idx==i)) = mycolors(i,2);
			I(:,:,2) = Itmp;
			Itmp = I(:,:,3);
			Itmp(myind(idx==i)) = mycolors(i,3);
			I(:,:,3) = Itmp;                      
	end

    s2 = 6;
    imshow(imenlarger(I, s2))
%     imshow(I)
    if savemyresults
       regime = I; %* 255;
       geotiffwrite([basedir 'clusters/regime_' num2str(i) '.tif'], regime, ref);
    end  
	if i <= nCluster
		%title(['Fire regime ' num2str(i)], 'FontSize', 20, 'FontName', 'Arial', 'FontWeight', 'Bold') % press any button
                fprintf('The number of grid cells in fire regime %d is %d\n', i, numel(myind(idx==i)));
		if savemyresults
		  print(hfig2, '-dpng', '-r500', [basedir 'clusters/regime-' num2str(i) '_' date '.png'])
		end            
		%pause
	end 
end

% take a look to the scores where you want to compute clusters
%
%insert other colors if you consider more than 4 clusters
%mycolors2 = {'m','r','b', 'g', 'c','y','k'}; %% the colors are showed exactly same as figure 135? I meant the corresponding color in each cluster (I guess so)
markers = {'o','s','+','x','d','*','v','h','^','p'};
type_of_imshow = 1;
use_color_image = 1;
switch type_of_imshow
  case 1 %2D
    figure(136)
    hold on
    for i = 1:nCluster
%       plot(XX(1,idx==i),XX(2,idx==i),[mycolors2{i} 'x'],'LineWidth', 2,'MarkerSize', 8)
%      plot(XX(1,idx==i),XX(2,idx==i),[mycolors2{i} markers{i}],'LineWidth', 2,'MarkerSize', 8)
        h=plot(XX(1,idx==i),XX(2,idx==i),[markers{i}],'LineWidth', 2,'MarkerSize', 8);
        switch use_color_image
		case 0
			I = zeros(m,n);
			I(myind(idx==i)) = 1;
		case 1
			I(:,:,1) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,2) = ones(m,n);%if you want black background use zeros intead of ones
			I(:,:,3) = ones(m,n);%if you want black background use zeros intead of ones
			
			Itmp = I(:,:,1);
			Itmp(myind(idx==i)) = mycolors(i,1);
			I(:,:,1) = Itmp;
			Itmp = I(:,:,2);
			Itmp(myind(idx==i)) = mycolors(i,2);
			I(:,:,2) = Itmp;
			Itmp = I(:,:,3);
			Itmp(myind(idx==i)) = mycolors(i,3);
			I(:,:,3) = Itmp;                      
        end
        if use_color_image == 1
            set(h,'color', mycolors(i,:));
        else
            set(h,'color','k');
        end
    end
    plot(regime_classes_centers(:,1),regime_classes_centers(:,2),'ko','LineWidth', 3,'MarkerSize', 12 )
    % multi-axes for 2D only
    if 0
      ah = copyobj(gca, gcf);
      set(ah, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'FontSize', 12, 'FontWeight', 'Bold');
    end
    %% title('Clusters based on principal components');
    legend({'Fire regime 1', 'Fire regime 2', 'Fire regime 3', 'Fire regime 4', 'Fire regime 5', 'Fire regime 6','Centroids'},'Location','northeast','FontSize',12)
    %legend({'Fire regime 1', 'Fire regime 2', 'Fire regime 3', 'Fire regime 4', 'Fire regime 5', 'Centroids'},'Location','northeast','FontSize',12) 
    %legend({'Fire regime 1', 'Fire regime 2', 'Fire regime 3', 'Fire regime 4', 'Centroids'},'Location','northeast','FontSize',12) 
    xlabel('Score of the 1st principal component', 'FontSize', 12);
    ylabel('Score of the 2nd principal component', 'FontSize', 12); %, 'FontWeight', 'Bold'
    set(gca, 'FontSize', 12, 'FontWeight', 'Bold' );
  case 2 % 3D
    figure(136) 
    for i = 1:nCluster
      scatter3(XX(1,idx==i), XX(2,idx==i), XX(3,idx==i),[markers{i}],'MarkerEdgeColor','k','MarkerFaceColor',mycolors(i,:));
      % disp(mycolors2{i});
      hold on;
    end
    xlabel('Score of the 1st principal component', 'FontSize', 12, 'FontWeight', 'Bold');
    ylabel('Score of the 2nd principal component', 'FontSize', 12, 'FontWeight', 'Bold');
    zlabel('Score of the 3rd principal component', 'FontSize', 12, 'FontWeight', 'Bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'Bold' );    
end

if savemyresults
	saveas(gcf,[basedir 'clusters/clusters_center_' date '.png'],'png')
end

%% One-way multivariate analysis of variance to test significance
if savemyresults
  for i=1:length(idx)
    group{i}=num2str(idx(i));
  end
  [d,p,stats]=manova1(XX', group);
end

