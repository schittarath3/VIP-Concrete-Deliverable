function newRepo = translateToPoint(aggRepo, stepSize, point)
%Translates all aggregates to defined point and checks for overlaps
%Inputs:
%   aggRepo: struct of aggregates with form of at least aggName->Points
%   stepSize: change of points for each step. Used 0.05.
%   point: point to which aggregates are translated
%Output:
%   newRepo: a struct with form of at least aggName->Points with translated
%            points
 
    %setting up constants 
    aggNames = fieldnames(aggRepo);
    numAggs = length(aggNames);
    
    %find distance of aggregate centroids to defined point
    newRepo = aggRepo;
    distTable = cell(numAggs, 2);
    for i = 1:numAggs
        curName = aggNames{i};
        aggCent = getCentroid(aggRepo.(curName).Points);
        aggDist = getDist(point, aggCent);
        distTable{i,1} = curName;
        distTable{i,2} = aggDist;
    end
    
    %sort array by ascending distance and get minimum and maximum
    distTable = sortrows(distTable, 2);
    minDist = min([distTable{:,2}]);
    maxDist = max([distTable{:,2}]);
    
    %normalize the closest aggregate to point which serve as the initial
    %overlap check
    newRepo.(distTable{1,1}).Points = normalizeTo(newRepo.(distTable{1,1}).Points, point);
    
    %iteratively translating aggregate to point
    for i = 2:numAggs
        state = false;
        curAggName = distTable{i,1};
        curAggPoints = newRepo.(curAggName).Points;
        
        %translate aggregate closer to point. scaleStep is normalized by
        %distance because greater distance are translated further than
        %closer aggregates
        while 1
            curCent = getCentroid(curAggPoints);
            curDist = getDist(point, curCent);
            %normalize travel distance
            if curDist < minDist
                curDist = minDist;
            end
            normScale = distNorm(maxDist, minDist, curDist);
            newCent = curCent + ((point - curCent) * stepSize * normScale);
            curAggPoints = normalizeTo(curAggPoints, newCent);
            
            %check with other aggregates using alphaShapes. should this be
            %the parfor loop?
            for j = 1:numAggs
                if j == i
                    continue
                end
                curOtherAgg = distTable{j,1};
                otherAggPoints = newRepo.(curOtherAgg).Points;
                otherAlpha = alphaShape(otherAggPoints);
                
                pointCheck = inShape(otherAlpha, curAggPoints);
                pointCheckSum = sum(pointCheck, 'all');
                if pointCheckSum > 0 %breaks if there is overlap and scales back once
                    curAggPoints = normalizeTo(curAggPoints, curCent);
                    newRepo.(curAggName).Points = curAggPoints(:,1:3);
                    state = true;
                    break
                end
            end
            if state == true
                break
            end
        end
    end
end

function normScale = distNorm(max, min, curDistance)
%normalize distance formula
%Inputs: 
%   max: maximum distance
%   min: minimum distance
%   curDiostance: current distance of aggregate to point
%Outputs:
%   normScale: double [0,1]
    normScale = (curDistance - min)/(max - min); 
    normScale = 1 - normScale;
end

function distance = getDist(pointFinal, pointsAggCent)
%gets distance to a point
%Inputs:
%   pointsFinal: final point
%   pointsAggCent: starting point, aggregate centroid
%Ouput:
%   distance
    distance = sqrt((pointFinal(1) - pointsAggCent(1))^2 + ...
                        (pointFinal(2) - pointsAggCent(2))^2 + ...
                          (pointFinal(3) - pointsAggCent(3))^2);
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

function datapointsn = normalizeTo(datapoints, newCentroid)
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
dcm = centroid - newCentroid; %Distance from centroid to cubeCentroid
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = [datapoints(vertice,:) - dcm];
end
end