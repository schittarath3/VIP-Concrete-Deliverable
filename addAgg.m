function newRepo = addAgg(originalRepo, aggRepo, cubeSize, newAggNum);
%Adds aggregates to the existing aggRepo
%Inputs: 
%   originalRepo: the repository generated using reducemesh.m
%   aggRepo: any repository generated after insertAgg.m
%   cubeSize: the max x, y, or z of the cube defined in cublets.m
%   newAggNum: the number of aggregates to be added to aggRepo
%Outputs:
%   newRepo: aggRepo with new aggregates added

    newRepo = aggRepo;
    
    %Setting up constants
    aggNames = fieldnames(originalRepo);
    numAgg = length(aggNames);
    oriNames = fieldnames(originalRepo.(aggNames{1}).Orientation);
    numOri = length(oriNames);
    
    %Adding aggregates
    for i = 1:newAggNum
        while 1
            aggInd = randperm(numAgg, 1); 
            oriInd = randperm(numOri, 1);

            %get a random aggregate and a orientation from originalRepo
            curAggName = aggNames{aggInd};
            curOriName = oriNames(oriInd);
            
            %get orientation index for storage
            oriX = str2num(curOriName{1}(8));
            oriY = str2num(curOriName{1}(16));
            oriZ = str2num(curOriName{1}(24));
            oriMat = [oriX oriY oriZ];
            
            %check to see if aggregate and orientation already exist in
            %aggRepo
            aggRepoNames = fieldnames(aggRepo);
            matches = 0;
            for arAggs = 1 : length(aggRepoNames)
                if strcmp(aggRepo.(aggRepoNames{i}).Original, curAggName) ... 
                        && sum(aggRepo.(aggRepoNames{i}).Orientation == oriMat) == 3
                    matches = matches + 1;
                end
                        
            end
            if matches == 0
                break
            elseif matches > 0
                continue
            end
        end
        
        %Select face of cube to start at
        switch (randperm(6, 1))
             case 1
                 newPoint = [randperm(cubeSize+1, 1)-1 0 randperm(cubeSize+1, 1)-1];
             case 2
                 newPoint = [cubeSize randperm(cubeSize+1, 1)-1 randperm(cubeSize+1, 1)-1];
             case 3
                 newPoint = [randperm(cubeSize+1, 1)-1 cubeSize randperm(cubeSize+1, 1)-1];
             case 4
                 newPoint = [0 randperm(cubeSize+1, 1)-1 randperm(cubeSize+1, 1)-1];
             case 5
                 newPoint = [randperm(cubeSize+1, 1)-1 randperm(cubeSize+1, 1)-1 cubeSize];
             case 6
                 newPoint = [randperm(cubeSize+1, 1)-1 randperm(cubeSize+1, 1)-1 0];
        end
        
        %Store new data into output repository
        curNewName = strcat("new", "_", curAggName, , curOriName{1}(8), ...
                            curOriName{1}(16), curOriName{1}(24));
        newRepo.(curNewName).Original = curAggName;
        newRepo.(curNewName).OriginalPoints = originalRepo.(curAggName).OriginalPoints;
        newRepo.(curNewName).OriginalFaces = originalRepo.(curAggName).OriginalFaces;
        newRepo.(curNewName).Orientation = oriMat;
        newRepo.(curNewName).Points = normalize(originalRepo.(curAggName).Vertices, newPoint);
        newRepo.(curNewName).Faces = originalRepo.(curAggName).Faces;
        
    end
end

function datapointsn = normalize(datapoints, newCentroid)
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