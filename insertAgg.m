function [insertRepo, aggRepo] = insertAgg(aggRepo, cubesCell, targetNum)
%Inserts aggergates into predefined cubes at the cube's centroid
%Inputs: 
%   aggRepo: Struct of aggregates in form (aggregateName-Faces,Orientations->orientations) 
%   cubesCell: Cell of cube vertice points
%   scaleFactor: A scaling factor for all aggregates recommended <0.01
%   numOrientations: The number of orientations that are allowed
%Ouput: 
%   insertRepo: Struct of aggregates in form
%   (aggregateName->Points,Faces,cubeNum). The number of aggregates is
%   number of aggregates * numOrientations

    %setting up constants
    numSlots = length(cubesCell); %number of available slots
    aggNames = fieldnames(aggRepo); %names of aggregates in AggRepo
    aggNames = aggNames(end-26:end);
    
    %creates new struct
    insertRepo = struct;
    cubeNum = randperm(27);
    if numSlots == targetNum
        for i = 1:targetNum
            insertRepo.(aggNames{i}) = aggRepo.(aggNames{i});
            insertRepo.(aggNames{i}).cubeNum = cubeNum(i);
        end
        aggRepo = rmfield(aggRepo, aggNames);
        
        newAggName = fieldnames(insertRepo); %get new aggregate names
        newRanInd = randperm(targetNum);
        for i = 1:targetNum %associate each aggregate in aggRepo to a cube in cubeCell
            curAggName = newAggName{newRanInd(i)}; 
            cubeNum = insertRepo.(curAggName).cubeNum;
            cubelet = cubesCell{cubeNum}; %gets random cube point
            cubeAlpha = alphaShape(cubelet);
            
            cubeCentroid = getCentroid(cubelet);
            insertRepo.(curAggName).Points ... %normalize aggregate centroid to cube centroid
                = normalizeTo(insertRepo.(curAggName).Points, cubeCentroid);
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

function datapointsn = normalizeTo(datapoints, cubeCentroid)
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
