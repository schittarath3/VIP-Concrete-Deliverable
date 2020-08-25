function sphereCell = VoidSearch(aggRepo, oriCubeSize, cubeSize, targetRate, spherePoints)
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
    
    %Growth rate
    scaleStepFactor = 0.8;
    
    %Getting the coverage rate needed as an end condition
    [covRate, totalVolume] = coverageRate(aggRepo, cubeSize^3);
    sphereCovRate = targetRate - covRate;
    
    %Creating the sphereCell
    sphereCell = cell(0,0);
    sphereNum = 0;
    curSphereRate = 0;
    
    %Looping until the sphere coverage rate is achieved
    while curSphereRate < sphereCovRate
        %Generates a new point
        xyz = getNewPos(oriCubeSize, cs);
        
        %Checking if point is within an aggregate
        if ~isPointIn(xyz, aggRepo)
            continue
        end
      
        %Create a sphere at the point
        sphere = normalizeTo(spherePoints, xyz);
        
        %Checking is sphere is out of bounds of the sub-cube
        if sum(inShape(cubeAlpha, sphere)) < length(sphere)
               continue
        end
         
        %Checking if the cube is already in another cube
        if ~inSphere(sphereCell, sphere)
            continue
        end

        growIter = 1;
        while 1
            %Grow sphere
            newSphere = growSphere(sphere, scaleStepFactor, xyz, growIter);
            %Performs a check if the sphere is within any sphere,
            %aggregate, or not within the sub-cube
            isInSomething = inCheck(cubeAlpha, sphereCell, aggRepo, newSphere);
            
            %Stops growing if the previous condition is true
            if ~isInSomething
                sphereNum = sphereNum + 1;
                %Saves to the output cell
                [sphereVolume, sphereCell] = saveToSphereCell(sphere, sphereCell, xyz, sphereNum, true, scaleStepFactor, growIter);
                curSphereRate = curSphereRate + (sphereVolume/cubeSize^3)
                break
            end
            growIter = growIter + 1
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

function [sphereVolume,sphereCell] = saveToSphereCell(sphere, sphereCell, xyz, sphereNum, grow, scaleStepFactor, growIter)
        
    if grow
            if growIter == 1
                scaleStepFactor = 1;
            else
                scaleStepFactor = scaleStepFactor * (growIter - 1);
            end
            scaleUp = [scaleStepFactor 0 0; 0 scaleStepFactor 0; 0 0 scaleStepFactor];
            sphere = normalizeTo(sphere, [0 0 0]);
            sphere = sphere * scaleUp;
            sphere = normalizeTo(sphere, xyz);
        end

        sphereLength = maxDiam(sphere);
        sphereVolume = (4/3)*pi*((sphereLength/2)^3);
        
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

function checksOut = inCheck(cubeAlpha, sphereCell, aggRepo, sphere)
%Checks if a sphere is within an aggregate or another sphere, and not
%within the sub-cube
    aggNames = fieldnames(aggRepo);
    aggNameLength = length(aggNames);
    checksOut = true;
    
    if sum(inShape(cubeAlpha, sphere)) < length(sphere)
       checksOut = false;
       return
    end
    if ~inSphere(sphereCell, sphere)
        checksOut = false;
        return
    end
    for agg = 1:aggNameLength
        curAggShape = alphaShape(aggRepo.(aggNames{agg}).Points);
        if  sum(inShape(curAggShape, sphere)) > 0
            checksOut = false;
            return
        end
    end
end

function checksOut = isPointIn(xyz, aggRepo)
%Checks if point is within a aggregate
        aggNames = fieldnames(aggRepo);
        aggNameLength = length(aggNames);
        checksOut = true;
        
        for agg = 1:aggNameLength
            curAggShape = alphaShape(aggRepo.(aggNames{agg}).Points);
            c = inShape(curAggShape, xyz(1), xyz(2), xyz(3));
            if inShape(curAggShape, xyz(1), xyz(2), xyz(3))
                checksOut = false;
                return
            end
        end
end

function c = growSphere(sphere, scaleStepFactor, xyz, growIter)
%Grows sphere by the specified scale factor
    scaleStepFactor = scaleStepFactor * growIter;
    scaleUp = [scaleStepFactor 0 0; 0 scaleStepFactor 0; 0 0 scaleStepFactor];
    a = normalizeTo(sphere, [0 0 0]);
    b = a * scaleUp;
    c = normalizeTo(b, xyz);
end

function checksOut = inSphere(sphereCell, sphere)
%Checks if sphere is within another sphere
        checksOut = true;
        if ~isempty(sphereCell)
            numSph = size(sphereCell);
            for sph = 1:numSph(1)
                spAlpha = alphaShape(sphereCell{sph});
                if sum(inShape(spAlpha, sphere)) > 0
                    checksOut = false;
                    return
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
        for row = 1:length(sphereCell)
            if round(sphereCell{row, 3}, 4) == bin
                num = num + 1;
                newSphereCell{i, 1}(num,:) = sphereCell(row,1:2);
            end
        end
    end
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