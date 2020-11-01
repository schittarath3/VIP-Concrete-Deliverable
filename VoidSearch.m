function sphereCell = VoidSearch(aggRepo, oriCubeSize, cubeSize, targetRate, grownSpheres)
%Attempts to find empty space in a given cube with or without aggregates. Randomly finds points that are within a sub-cube of the original cube and grows spheres at the empty point
%Inputs:
%   aggRepo: a struct of form of at least (aggName -> Points), can be empty
%   oriCubeSize: the original cube size in which the aggregates
%   cubeSize: the sub-cube where empty points will be located
%   targetRate: the target coverage rate
%   spherePoints: the points of an stl file that will be used to imitate an
%   aggregate as a cell
%Outputs:
%   sphereCell: a cell of form (sphereLength -> Points, centroid)
%   containing all spheres generated
    
    %Creating the sub-cube alphashape to check for inShape
    cs = cubeSize;
    cubePoints = [0 0 0; 0 cs 0; cs 0 0; ...
    cs cs 0; 0 0 cs; 0 cs cs; ...
    cs 0 cs; cs cs cs];
    cubePoints = normalizeTo(cubePoints, [oriCubeSize/2 oriCubeSize/2 oriCubeSize/2]);
    cubeAlpha = alphaShape(cubePoints);
    
    %Getting the coverage rate needed as an end condition
    [covRate, totalVolume] = coverageRate(aggRepo, cubeSize^3);
    sphereCovRate = targetRate - covRate
    
    %Creating the sphereCell
    sphereCell = cell(0,0);
    sphereNum = 0;
    curSphereRate = 0;
    tryNum = 0;
    
%Looping until the sphere coverage rate is achieved
while curSphereRate < sphereCovRate
    if tryNum > 1000
        disp("volume fraction cannot be reached")
        break
    end
    tryNum = tryNum + 1;
    %Generates a new point
    xyz = getNewPos(oriCubeSize, cs);
    

    %Checking if point is within an aggregate
    if ~isPointIn(xyz, aggRepo)
        continue
    end
    
    %Create a sphere at the point
    sphere = normalizeTo(grownSpheres{1}, xyz);

    %Checking is sphere is out of bounds of the sub-cube
    if sum(inShape(cubeAlpha, sphere)) < length(sphere)
           continue
    end

    %Checking if the sphere is already in another sphere
    if ~inSphere(sphereCell, sphere)
        continue
    end

    growIter = 2;
    while 1
        tryNum = 0;
        if growIter == 4
            sphereNum = sphereNum + 1;
            [sphereVolume, sphereCell] = saveToSphereCell(newSphere, sphereCell, xyz, sphereNum);
            break;
        end
        %Grow sphere
        newSphere = normalizeTo(grownSpheres{growIter}, xyz);
        %Performs a check if the sphere is within any sphere,
        %aggregate, or not within the sub-cube
        isInSomething = inCheck(cubeAlpha, sphereCell, aggRepo, newSphere, xyz);

        %Stops growing if the previous condition is true
        if ~isInSomething
            sphereNum = sphereNum + 1;
            %Saves to the output cell
            [sphereVolume, sphereCell] = saveToSphereCell(grownSpheres{growIter-1}, sphereCell, xyz, sphereNum);
            curSphereRate = curSphereRate + (sphereVolume/cubeSize^3);
            break
        end
        growIter = growIter + 1;
    end
end
    %Sort the sphere cell by length
    sphereCell = sortSphereCell(sphereCell);
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

function [sphereVolume,sphereCell] = saveToSphereCell(sphere, sphereCell, xyz, sphereNum)
%     if growIter > 2
%         scaleStepFactor = 1 + scaleStepFactor * (growIter - 1);
%         sphere = normalizeTo(sphere, [0 0 0]);
%         sphere = Scale(sphere, scaleStepFactor, scaleStepFactor, scaleStepFactor);
%         sphere = normalizeTo(sphere, xyz);
%     end

    sphere = normalizeTo(sphere, xyz);
    sphereLength = maxDiam(sphere);
    sphereVolume = volume(alphaShape(sphere));

    sphereCell{sphereNum, 1} = sphere;
    sphereCell{sphereNum, 2} = xyz;
    sphereCell{sphereNum, 3} = sphereLength;
end

function xyz = getNewPos(oriCubeSize, cs)
%Generates a random point within the sub-cube
    maxXYZ = oriCubeSize/2 + round(cs*0.90);
    minXYZ = oriCubeSize/2 - round(cs*0.90);
    xPos = round((maxXYZ -  minXYZ).*rand +  minXYZ);
    yPos = round((maxXYZ -  minXYZ).*rand +  minXYZ);
    zPos = round((maxXYZ -  minXYZ).*rand +  minXYZ);
    xyz = [xPos yPos zPos];
end

function checksOut = inCheck(cubeAlpha, sphereCell, aggRepo, sphere, xyz)
%Checks if a sphere is within an aggregate or another sphere, and not
%within the sub-cube
    aggNames = getClosestAggregates(xyz, aggRepo);
    aggNameLength = length(aggNames);
    checksOut = true;
    
    if sum(inShape(cubeAlpha, sphere)) < length(sphere)
       checksOut = false;
       return
    end
    
    if ~isempty(sphereCell) && ~inSphere(sphereCell, sphere)
        checksOut = false;
        return
    end
    
    for agg = 1:aggNameLength
        crit = criticalAlpha(alphaShape(aggRepo.(aggNames{agg}).Points), 'one-region') + 20;
        curAggShape = alphaShape(aggRepo.(aggNames{agg}).Points, crit);
        if  sum(inShape(curAggShape, sphere)) > 0
            checksOut = false;
            return
        end
    end
    
end

function checksOut = isPointIn(xyz, aggRepo)
%Checks if point is within a aggregate
    aggNames = getClosestAggregates(xyz, aggRepo);
    aggNameLength = length(aggNames);
    checksOut = true;

    for agg = 1:aggNameLength
        curAggShape = alphaShape(aggRepo.(aggNames{agg}).Points);
        if inShape(curAggShape, xyz(1), xyz(2), xyz(3))
            checksOut = false;
            return
        end
    end
end

function c = growSphere(sphere, scaleStepFactor, xyz, growIter)
%Grows sphere by the specified scale factor
    scaleStepFactor = 1 + (scaleStepFactor * growIter);
    a = normalizeTo(sphere, [0 0 0]);
    b = Scale(a, scaleStepFactor, scaleStepFactor, scaleStepFactor);
    c = normalizeTo(b, xyz);
end

function checksOut = inSphere(sphereCell, sphere)
%Checks if sphere is within another sphere
    checksOut = true;
    if ~isempty(sphereCell)
        sphC = getClosestSpheres(getCentroid(sphere), sphereCell);
        numSph = size(sphC);
        for sph = 1:numSph(1)
            crit = criticalAlpha(alphaShape(sphC{sph}), 'one-region') + 20;
            spAlpha = alphaShape(sphC{sph}, crit);
            if sum(inShape(spAlpha, sphere)) ~= 0
                checksOut = false;
                break
            elseif sum(inShape(alphaShape(sphere), sphC{sph})) ~= 0
                checksOut = false;
                break
            end
        end
    end
end

function newSphereCell = sortSphereCell(sphereCell)
%Sorts all the spheres by length
    bins = unique(round(cell2mat(sphereCell(:,3)), 4));
    numBins = length(bins);
    
    newSphereCell = cell(numBins,2);
    newSphereCell(:,2) = num2cell(bins);
    for i = 1:numBins
        bin = bins(i);
        newSphereCell{i, 1} = cell(0,2);
        num = 0;
        rowNum = size(sphereCell);
        for row = 1:rowNum(1)
            if round(sphereCell{row, 3}, 4) == bin
                num = num + 1;
                newSphereCell{i, 1}(num,:) = sphereCell(row,1:2);
            end
        end
    end
end

function closestSpheres = getClosestSpheres(xyz, sphereCell)
    if isempty(sphereCell)
        return
    end
    sphereNum = size(sphereCell);
    sphereNum = sphereNum(1);
    tempCell = sphereCell;
    for i = 1:sphereNum
        tempCell(i,4) = mat2cell(tempCell{i,2} - xyz, [1]);
    end
    tempCell(:,4) = num2cell(cellfun(@norm, tempCell(:,4)));
    tempCell = sortrows(tempCell, 4);
    
    if sphereNum < 15
         closestSpheres = tempCell(1:sphereNum, 1);
    else
        closestSpheres = tempCell(2:15, 1);
    end
end

function closestAggs = getClosestAggregates(xyz, aggRepo)
    aggs = fieldnames(aggRepo);
    for i = 1:length(aggs)
        aggs(i,2) = mat2cell(norm(getCentroid(aggRepo.(aggs{i}).Points) - xyz), [1]);
    end
    aggs = sortrows(aggs, 2);
    closestAggs = aggs(1:10);
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

function dist = maxDiam(aggpts)
dist = 0;
for i = 1:length(aggpts)
    for k = 1:length(aggpts)
    distp = abs(norm(aggpts(i,:) - aggpts(k,:)));
    if distp > dist
        dist = distp;
    end
    end
end
end

function nv = Scale(pts,sx,sy,sz) 
%Rotational matrix
scalem = diag([sx sy sz]);

nv = pts * scalem;
end