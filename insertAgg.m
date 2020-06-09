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
    
    insertRepo = aggRepo;
    if numSlots == totalNumAgg
        for i = 1:numAggs %Rescaling points by scaleFactor
            oriRandIndex = randperm(totalOrientations, numOrientations);
            for x = 1:numOrientations
                oriName = orientationNames{oriRandIndex(x)}
                insertRepo.(aggNames{i}).newOrientations.(oriName)...
                    = insertRepo.(aggNames{i}).Orientation.(oriName) * scaleMat;
            end
        end
        for i = 1:numAggs %associate each aggregate in aggRepo to a cube in cubeCell
            aggName = aggNames{aggRandIndex(i)}
            newOriNames = fieldnames(insertRepo.(aggName).newOrientations);
            for x = 1:numOrientations
                cubeNum = cubeRandIndex(i);
                cubelet = cubesCell{cubeNum};
                cubeAlpha = alphaShape(cubelet);
                insertRepo.(aggName).cubeNum = cubeNum;
                cubeCentroid = getCentroid(cubelet);
                curOriName = newOriNames{x}
                insertRepo.(aggName).newOrientations.(curOriName)... 
                    = normalize(insertRepo.(aggName).newOrientations.(curOriName), cubeCentroid);
                pointCheck = ~inShape(cubeAlpha, insertRepo.(aggName).newOrientations.(curOriName));
                if pointCheck > 0
                    disp("Aggregate " + aggName + curOriName +...
                         " does not fit within mini cube. Try decreasing the scaling factor.")
                    break
                else
                    continue
                end
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
