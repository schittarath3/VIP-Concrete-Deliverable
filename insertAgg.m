function insertRepo = insertAgg(aggRepo, cubesCell, targetNum)
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
    numAggs = length(aggNames); %number of aggregates
    
    %more constants
    orientationNames = fieldnames(aggRepo.(aggNames{1}).Orientation); %names of orientations
    totalOrientations = length(orientationNames);
    
    %creates new struct
    insertRepo = struct;
    if numSlots == targetNum
        curNumAgg = 1;
        cubeInd = 1;
        while curNumAgg <= targetNum
            aggRandNum = randi(targetNum, 1);
            curAggName = aggNames{aggRandNum};
            oriRandNum = randi(totalOrientations, 1);
            oriName = orientationNames{oriRandNum};
                
            %Create index for orientation
            oriX = str2num(oriName(8));
            oriY = str2num(oriName(16));
            oriZ = str2num(oriName(24));
            oriMat = [oriX oriY oriZ];
            
            newName = strcat(curAggName, "_", oriName(8), oriName(16), ...
                                oriName(24));
                            
            if isfield(insertRepo, newName)
                continue
            end
                            
           insertRepo.(newName).Original = aggNames{aggRandNum}; %Store original number
           insertRepo.(newName).cubeNum = cubeInd; %associate aggregate with cubeNum
           insertRepo.(newName).OriginalPoints = aggRepo.(aggNames{aggRandNum}).OriginalPoints;
           insertRepo.(newName).OriginalFaces = aggRepo.(aggNames{aggRandNum}).OriginalFaces;
           cubeInd = cubeInd + 1;
           insertRepo.(newName).Orientation = oriMat; %Indices for linspace
           insertRepo.(newName).Points ... 
                = aggRepo.(aggNames{aggRandNum}).Orientation.(oriName);
           insertRepo.(newName).Faces = aggRepo.(aggNames{aggRandNum}).Faces; %Store connectivity
           
           curNumAgg = curNumAgg + 1;
        end

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
            pointCheck = ~inShape(cubeAlpha, insertRepo.(curAggName).Points);
            pointCheckSum = sum(pointCheck, 'all');
            if pointCheckSum > 0 %checks if any aggregate points are outside the cube
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
