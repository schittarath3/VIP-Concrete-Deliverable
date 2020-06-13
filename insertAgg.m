function insertRepo = insertAgg(aggRepo, cubesCell, scaleFactor, numOrientations)

    %setting up constants
    numSlots = length(cubesCell);
    numAggs = length(fieldnames(aggRepo));
    totalNumAgg = numAggs * numOrientations;
    scaleMat = [scaleFactor 0 0; 0 scaleFactor 0; 0 0 scaleFactor];
    
    cubeRandIndex = randperm(numSlots); %psuedorandom array to index cubes
    aggRandIndex = randperm(numAggs); %psuedorandom array to index aggregates
    aggNames = fieldnames(aggRepo);
    orientationNames = fieldnames(aggRepo.(aggNames{1}).Orientation);
    totalOrientations = length(orientationNames);
    
    insertRepo = struct;
    if numSlots == totalNumAgg
        for i = 1:numAggs %Rescaling points by scaleFactor
            aggRandNum = aggRandIndex(i);
            oriRandIndex = randperm(totalOrientations, numOrientations);
            for x = 1:numOrientations
                oriName = orientationNames{oriRandIndex(x)};
                newName = strcat(aggNames{aggRandNum}, oriName);
                insertRepo.(newName).Points ...
                    = aggRepo.(aggNames{aggRandNum}).Orientation.(oriName) * scaleMat;
                insertRepo.(newName).Faces = aggRepo.(aggNames{aggRandNum}).Faces;
            end
        end
        
        newAggName = fieldnames(insertRepo);
        for i = 1:totalNumAgg %associate each aggregate in aggRepo to a cube in cubeCell
            curAggName = newAggName{i};
            cubeNum = cubeRandIndex(i);
            cubelet = cubesCell{cubeNum};
            cubeAlpha = alphaShape(cubelet);
            insertRepo.(curAggName).cubeNum = cubeNum;
            cubeCentroid = getCentroid(cubelet);
            insertRepo.(curAggName).Points ... 
                = normalize(insertRepo.(curAggName).Points, cubeCentroid);
            pointCheck = ~inShape(cubeAlpha, insertRepo.(curAggName).Points);
            pointCheckSum = sum(pointCheck, 'all');
            if pointCheckSum > 0
                disp("Aggregate " + curAggName +...
                     " does not fit within mini cube. Try decreasing the scaling factor.")
                break
            else
                continue
            end
        end
    else
        disp("Number of cubes and number of aggregates do not match")
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

function datapointsn = normalize(datapoints, cubeCentroid)
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
dcm = centroid - cubeCentroid; %Distance from centroid to cubeCentroid
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = [datapoints(vertice,:) - dcm];
end
end
