function insertRepo = insertAgg(aggRepo, cubesCell, scaleFactor)

%setting up constants
    numSlots = length(cubesCell);
    numAggs = length(fieldnames(aggRepo));
    scaleMat = [scaleFactor 0 0; 0 scaleFactor 0; 0 0 scaleFactor];
    cubeRandIndex = randperm(numSlots); %psuedorandom array to index cubes
    aggRandIndex = randperm(numAggs); %psuedorandom array to index aggregates
    
    aggNames = fieldnames(aggRepo);
    insertRepo = aggRepo;
    if numSlots == numAggs
        for i = 1:numAggs %Rescaling points by scaleFactor
            insertRepo.(aggNames{i}).Points = insertRepo.(aggNames{i}).Points * scaleMat;
        end
        for i = 1:numAggs %associate each aggregate in aggRepo to a cube in cubeCell
            cubeNum = cubeRandIndex(i);
            cubelet = cubesCell{cubeNum};
            cubeAlpha = alphaShape(cubelet);
            aggName = aggNames{aggRandIndex(i)};
            insertRepo.(aggName).cubeNum = cubeNum;
            cubeCentroid = getCentroid(cubelet);
            insertRepo.(aggName).Points = normalize(insertRepo.(aggName).Points, cubeCentroid);
            pointCheck = ~inShape(cubeAlpha, insertRepo.(aggName).Points);
            if pointCheck > 0
                disp("Aggregate " + aggName + " does not fit within mini cube. Try decreasing the scaling factor.")
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
