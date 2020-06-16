function insertRepo = insertAgg(aggRepo, cubesCell, scaleFactor, numOrientations)
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
    numAggs = length(fieldnames(aggRepo)); %number of aggregates
    totalNumAgg = numAggs * numOrientations; %total number of possible new aggregates
    scaleMat = [scaleFactor 0 0; 0 scaleFactor 0; 0 0 scaleFactor]; %scaling matrix
    
    %more constants
    aggRandIndex = randperm(numAggs); %psuedorandom array to index aggregates
    aggNames = fieldnames(aggRepo); %names of aggregates in AggRepo
    orientationNames = fieldnames(aggRepo.(aggNames{1}).Orientation); %names of orientations
    totalOrientations = length(orientationNames);
    
    %creates new struct
    insertRepo = struct;
    if numSlots == totalNumAgg
        for i = 1:numAggs %Rescaling points by scaleFactor
            aggRandNum = aggRandIndex(i); %getting random aggregate
            oriRandIndex = randperm(totalOrientations, numOrientations); %fetching random orientation
            
            %stores aggregates with orientation name
            %scaling aggregates
            for x = 1:numOrientations
                oriName = orientationNames{oriRandIndex(x)};
                newName = strcat(aggNames{aggRandNum}, num2str(i));
                
                %Create index for orientation
                oriX = str2num(oriName(8));
                oriY = str2num(oriName(16));
                oriZ = str2num(oriName(24));
                oriMat = [oriX oriY oriZ];
                
                insertRepo.(newName).Orientation = oriMat;
                insertRepo.(newName).Points ... 
                    = aggRepo.(aggNames{aggRandNum}).Orientation.(oriName) * scaleMat;
                insertRepo.(newName).Faces = aggRepo.(aggNames{aggRandNum}).Faces;
            
            end
        end
        
        newAggName = fieldnames(insertRepo); %get new aggregate names
        
        for i = 1:totalNumAgg %associate each aggregate in aggRepo to a cube in cubeCell
            curAggName = newAggName{i}; 
            cubeNum = i;
            cubelet = cubesCell{cubeNum}; %gets random cube point
            cubeAlpha = alphaShape(cubelet);
            
            insertRepo.(curAggName).cubeNum = cubeNum; %associate aggregate with cubeNum
            cubeCentroid = getCentroid(cubelet);
            insertRepo.(curAggName).Points ... %normalize aggregate centroid to cube centroid
                = normalize(insertRepo.(curAggName).Points, cubeCentroid);
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
