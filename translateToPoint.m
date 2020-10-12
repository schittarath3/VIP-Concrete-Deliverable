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
    minRange = (point - 5);
    maxRange = (point + 5);
    
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
    
    %normalize the closest aggregate to point which serve as the initial
    %overlap check
    newRepo.(distTable{1,1}).Points = normalizeTo(newRepo.(distTable{1,1}).Points, point);
    
    %iteratively translating aggregate to point
    for i = 2:numAggs
        state = false;
        curAggName = distTable{i,1};
        curAggPoints = newRepo.(curAggName).Points;
        step = 1;
        curCent = getCentroid(curAggPoints);
        curDist = getDist(point, curCent);
        if curDist < stepSize
            curDist = 10;
        end
        normScale = distNorm(curDist, stepSize);

        %translate aggregate closer to point. scaleStep is normalized by
        %distance because greater distance are translated further than
        %closer aggregates
        while 1
            %normalize travel distance
            newCent = curCent + (point - curCent)*normScale*step;
            curAggPoints = normalizeTo(curAggPoints, newCent);
            %check with other aggregates using alphaShapes. 
            closestAggs = getClosestAggregates(curCent, newRepo)';
            for j = 1:length(closestAggs)
                curOtherAgg = closestAggs{j,1};
                if strcmp(curAggName, curOtherAgg)
                    continue
                end
                otherAggPoints = newRepo.(curOtherAgg).Points;
                otherAlpha = alphaShape(otherAggPoints);
                crit = criticalAlpha(otherAlpha, 'one-region') + 40;
                otherAlpha = alphaShape(otherAggPoints, crit);
                
                pointCheck = inShape(otherAlpha, curAggPoints);
                pointCheckSum = sum(pointCheck, 'all');
                if pointCheckSum > 0 %breaks if there is overlap and scales back once
                    state = true;
                    break
                end
            end
            if sum(minRange < newCent) == 3 && sum(maxRange > newCent) == 3
                state = true;
            end
            if state == true
                if step ~= 1
                    newCent = curCent + (point - curCent)*normScale*(step - 1);
                    curAggPoints = normalizeTo(curAggPoints, newCent);
                    newRepo.(curAggName).Points = curAggPoints(:,1:3);
                    curAggOriPoints = newRepo.(curAggName).OriginalPoints;
                    curAggOriPoints = normalizeTo(curAggOriPoints, newCent);
                    newRepo.(curAggName).OriginalPoints = curAggOriPoints(:,1:3);
                end
                break
            end
            step = step + 1;
        end
    end
end

function normScale = distNorm(curDistance, stepSz)
%normalize distance formula
%Inputs: 
%   max: maximum distance
%   min: minimum distance
%   curDiostance: current distance of aggregate to point
%Outputs:
%   normScale: double [0,1]
    curDistStep = curDistance/stepSz;
    normScale = 1/round(curDistStep);
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

function closestAggs = getClosestAggregates(xyz, aggRepo)
    aggs = fieldnames(aggRepo);
    for i = 1:length(aggs)
        aggs(i,2) = mat2cell(norm(getCentroid(aggRepo.(aggs{i}).Points) - xyz), [1]);
    end
    aggs = sortrows(aggs, 2);
    closestAggs = aggs(1:30);
end 