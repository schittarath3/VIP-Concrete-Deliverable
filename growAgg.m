function newRepo = growAgg(aggRepo, cubeCell, scaleStepFactor)
%Grows aggregates within its respective cubes
%Input: Aggregate repo (aggRepo), Cell of cube coordinate (cubeCell), and a
%scaling factor >1.000 (scaleStepFactor). Ideally a small step size, >1.000
%<1.05. A smaller value means longer computational time
%Output: A new repoistory (struct) containing aggregate names and their
%repsecitve cube number and points

    %Setting up constants
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    
    %3D scaling matrices
    scaleUp= [scaleStepFactor 0 0 ;0 scaleStepFactor 0; 0 0 scaleStepFactor];
    scaleDown =[1/scaleStepFactor 0 0 ;0 1/scaleStepFactor 0; 0 0 1/scaleStepFactor];
    
    newRepo = struct %create new struct to store new points
    for i = 1:numAgg %iterates through aggRepo 
        
        %setting up variables
        curAggName = aggNames{i};
        curCubeNum = aggRepo.(curAggName).cubeNum;
        
        newRepo.(curAggName).Orientation = aggRepo.(curAggName).Orientation; %Store orientation
        newRepo.(curAggName).cubeNum = curCubeNum; %store cubeNum to new repo
        newRepo.(curAggName).Faces = aggRepo.(curAggName).Faces; %store connectivity
        
        curCubePoints = cubeCell{curCubeNum}; %get aggregate's points
        curCubeAlpha = alphaShape(curCubePoints); %get cube alphaShape
        curCubeCent = getCentroid(curCubePoints); %get cube centroid
        
        curPoints = aggRepo.(curAggName).Points; %setting up points for iteration
       
        while 1
            curPoints = normalize(curPoints, [0 0 0]); %translate points to origin
                                                       %for scaling
            curPoints = curPoints * scaleUp; %growing aggregate
            curPoints = normalize(curPoints, curCubeCent); %translate back to 
                                                           %cube centroid
            
            %checking for overlap with cube alphaShape (ie points are within)                                               
            pointCheck = ~inShape(curCubeAlpha, curPoints);
            pointCheckSum = sum(pointCheck, 'all');
            if pointCheckSum > 0 %breaks if there is overlap and scales back once
                curPoints = normalize(curPoints, [0 0 0]);
                curPoints = curPoints * scaleDown;
                curPoints = normalize(curPoints, curCubeCent);
                newRepo.(curAggName).Points = curPoints(:,1:3);
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

