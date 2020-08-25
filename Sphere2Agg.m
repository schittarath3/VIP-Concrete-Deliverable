function [newaggRepo, agg_dist]= Sphere2Agg(aggRepo,sphereCell,results,sieveSz,totalAggs)
%Pack the rest of the aggregates to the original 27 following the grain
%size distribution obtained from the 2D image analysis:
%Inputs:
%   aggRepo - the generated distributed repository containing the scaled
%   aggregates to fit inside the container.
%   sphereCell - the generated repository containing possible empty spaces
%   to fit the aggregates inside of the original packed 27 aggregates.
%   results - a mx1 cell containg the lengths (bins) of the grains 
%   sieveSz - a 1xn cell containing the results of the grain size
%   totalAggs - total number of aggregates contained

%Determine how many aggregates can be packed to get the desired
%distribution based on total number of aggs (user-determined):
agg_dist = floor(totalAggs*(results/100));

%Determine the total number of available sphere cells to fill
totasph_bins = size(sphereCell,1);
available_sph = zeros(1,totasph_bins);
for r = 1:totasph_bins
    available_sph(1,r) = size(sphereCell{r,1},1);
end

%Obtaining the threshold sieve size for each of the aggregates
binSz = cell2mat(sphereCell(:,2)); 
maxSz = binSz(end);
maxBin = find(sieveSz <= maxSz);
maxBin = maxBin(end);

%Sorting the aggregates repository by bin size and diameter
fields_agg = fieldnames(aggRepo);
for agg = 1:length(fields_agg)
    diameter = aggRepo.(fields_agg{agg}).Diameter;
    bin = aggRepo.(fields_agg{agg}).bin;
    aggRepoSort{bin,1}.(fields_agg{agg}).Index = agg;
    aggRepoSort{bin,1}.(fields_agg{agg}).Diameter = diameter;
    
    %Finding the sphere cell bin that matches with the aggregates
    fit = find(binSz>=diameter);
        if ~isempty(fit)
        sphereCell_fit = fit(1);
        aggRepoSort{bin,1}.(fields_agg{agg}).Spherebin = sphereCell_fit;
        end
end

for bins = 1:maxBin
    try
    for aggs = 1:agg_dist(bins)
        binSz_fields = fieldnames(aggRepoSort{bins,1});
        agg_idx = aggRepoSort{bins,1}.(char(binSz_fields(aggs))).Index;
        agg_sph_bin = aggRepoSort{bins,1}.(char(binSz_fields(aggs))).Spherebin;
        
        pts = aggRepo.(fields_agg{agg_idx}).OriginalPoints;
        sph_cm = cell2mat(sphereCell{agg_sph_bin,1}(1,2));

        %Removing the sphere cell used
        sphereCell{bins,1}(1:end-1,1:2) = sphereCell{bins,1}(2:end,1:2);
        sphereCell{bins,1}(end,1:2) = cell(1,2);

        %Rewriting the coordinates for the aggregates to substitute...
         newaggRepo.(fields_agg{agg_idx}).Original = aggRepo.(fields_agg{agg_idx}).Original;
         newaggRepo.(fields_agg{agg_idx}).OriginalPoints = translate2Center(pts, sph_cm);
         newaggRepo.(fields_agg{agg_idx}).OriginalFaces = aggRepo.(fields_agg{agg_idx}).OriginalFaces;
         newaggRepo.(fields_agg{agg_idx}).Diameter = aggRepo.(fields_agg{agg_idx}).Diameter;
    end
    
    catch
        disp(['Fail to fit bin ' num2str(bins) '--Consider lowering the number of aggregates'])
    end
end
end

function datapointsn = translate2Center(datapoints,center)
%Obtain the center of each aggregate
x = datapoints(:,1);
y = datapoints(:,2);
z = datapoints(:,3);

xcm = sum(x)./length(x);
ycm = sum(y)./length(y);
zcm = sum(z)./length(z);
centroid = [xcm ycm zcm];

%Obtain the matrix with the distance of each vertices to the center
datapointsn = datapoints;
dcm = center - centroid; %Distance from centroid to center
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = datapoints(vertice,:) + dcm;
end
end