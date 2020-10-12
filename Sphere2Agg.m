function [newaggRepo, agg_dist, aggRepo]= Sphere2Agg(aggRepo,sphereCell,results,sieveSz,iter,totalAggs)
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
%Outputs:
    %   newaggRepo - new repository of the aggregates containing the added
    %   scaled aggregates inside sphere cells
    %   agg_dist - distribution matrix containing the number of aggregates
    %   from each bin/sieve size packed
    %   sphereCell - updated sphere cell with the remaining sphere cells 

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

if length(aggRepoSort) < maxBin
    maxBin = length(aggRepoSort);
end

numRemove = 1;
removedFields = {};
for bins = 1:maxBin
    if isempty(aggRepoSort{bins,1})
        continue
    end
    binSz_fields = fieldnames(aggRepoSort{bins,1});
    if length(binSz_fields) < agg_dist(bins)
        numAgg = length(binSz_fields);
    else
        numAgg = agg_dist(bins);
    end
    for aggs = 1:numAgg
        agg_idx = aggRepoSort{bins,1}.(char(binSz_fields(aggs))).Index;
        agg_sph_bin = aggRepoSort{bins,1}.(char(binSz_fields(aggs))).Spherebin;
        
        pts = aggRepo.(fields_agg{agg_idx}).OriginalPoints;
        try
            sph_cm = cell2mat(sphereCell{agg_sph_bin,1}(1,2));
        catch
            continue
        end

        %Removing the sphere cell used
        sphereCell{agg_sph_bin,1} = sphereCell{agg_sph_bin,1}(2:end,1:2);
        newName = strcat(fields_agg{agg_idx}, "_", num2str(iter));
        %Rewriting the coordinates for the aggregates to substitute...
         newaggRepo.(newName).Original = aggRepo.(fields_agg{agg_idx}).Original;
         newaggRepo.(newName).OriginalPoints = translate2Center(pts, sph_cm);
         newaggRepo.(newName).OriginalFaces = aggRepo.(fields_agg{agg_idx}).OriginalFaces;
         newaggRepo.(newName).Points = translate2Center(aggRepo.(fields_agg{agg_idx}).Points, sph_cm);
         newaggRepo.(newName).Faces = aggRepo.(fields_agg{agg_idx}).Faces;
         newaggRepo.(newName).Diameter = aggRepo.(fields_agg{agg_idx}).Diameter;
         newaggRepo.(newName).Centroid = sph_cm;
         removedFields{numRemove, 1} = fields_agg{agg_idx};
         numRemove = numRemove + 1;
    end
end
aggRepo = rmfield(aggRepo, removedFields);
end

function centroid = getCentroid(datapoints)
    x = datapoints(:,1);
    y = datapoints(:,2);
    z = datapoints(:,3);

    xcm = sum(x)./length(x);
    ycm = sum(y)./length(y);
    zcm = sum(z)./length(z);
    centroid = [xcm ycm zcm];
end

function datapointsn = translate2Center(datapoints,center)
%Translate the aggregates to the origin (0,0,0)
%Inputs: 
    %datapoints - coordinates of the aggregates with each column
    %representing x,y,z
    %center - the point to translate the CoM of the aggregates towards
%Outputs:   
    %datapointsn - translated coordinates
    
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