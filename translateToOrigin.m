function [newRepo, distTable] = translateToOrigin(aggRepo, stepSize)
    
    aggNames = fieldnames(aggRepo);
    numAggs = length(aggNames);
    
    newRepo = aggRepo;
    distTable = cell(numAggs, 2);
    for i = 1:numAggs
        curName = aggNames{i};
        aggCent = getCentroid(aggRepo.(curName).Points);
        aggDist = sqrt(aggCent(1)^2 + aggCent(2)^2 + aggCent(3)^2);
        distTable{i,1} = curName;
        distTable{i,2} = aggDist;
    end
    
    distTable = sortrows(distTable, 2);
    
    newRepo.(distTable{1,1}).Points = normalizeTo(newRepo.(distTable{1,1}).Points, [0 0 0]);
    for i = 2:numAggs
        state = false;
        curAggName = distTable{i,1};
        curAggPoints = newRepo.(curAggName).Points;
        while 1
            curCent = getCentroid(curAggPoints);
            newCent = curCent + (([0 0 0] - curCent) * 0.05);
            curAggPoints = normalizeTo(curAggPoints, newCent);
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