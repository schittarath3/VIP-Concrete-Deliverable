clc
clear
%% pull files from directory
%pull container info and normalize
testContainer = stlread('STL Files\container25mm.stl');
testContainerNorm = normalize(testContainer.Points);
containerAlpha = alphaShape(testContainer.Points)

%pull aggregate files from directory
fileList = dir('STL Files\Aggregate Test\*.stl');
fileList = struct2cell(fileList);
fileListSize = size(fileList);
aggList = cell(1, fileListSize(2));
samples = size(fileList);
meshvol = zeros(1,samples(2));

for i = 1:samples(2)
    %obtain stl
    folder = fileList{2, i};
    file = fileList{1, i};
    aggName = file(1:end-4);
    fn = strcat(folder, '\', file);
    aggList{i} = stlread(fn);
   
    %generating repository of aggregates
    repos.(aggName).Points = normalize(aggList{i}.Points);
    repos.(aggName).ConnectivityList = aggList{i}.ConnectivityList;
end
%% testing stuff
numSamples = size(fileList);
numSamples = numSamples(1,2);
initAgg = fileList(1,randi(numSamples));
%initAggName = initAgg{1}(1:end-4)
initAggName = 'ellipsoid'


%setting up boundaries of container
% containerMaxX = max(testContainer.Points(:,1));
% containerMaxY = max(testContainer.Points(:,2));
% containerMaxZ = max(testContainer.Points(:,3));
% containerMinX = min(testContainer.Points(:,1));
% containerMinY = min(testContainer.Points(:,2));
containerMinZ = min(testContainer.Points(:,3))


%make min z-coord equal to container's min z-coord
zMinDist = min(repos.(initAggName).Points(:,3)) - containerMinZ - 1;
for i = 1:length(repos.(initAggName).Points(:,3))
    repos.(initAggName).Points(i,3) = repos.(initAggName).Points(i,3) - zMinDist;
end

%random translation of coordinates on xy plane
xTheta = (-1 + 2.*rand(1,1))*pi;
yTheta = (-1 + 2.*rand(1,1))*pi;
beforeTranslation = repos.(initAggName).Points;
xyTranslate(repos.(initAggName).Points, containerAlpha, xTheta, yTheta);
afterTranslation = repos.(initAggName).Points;



%% functions

function datapointsn = normalize(datapoints)
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
dcm = [0, 0, 0] - centroid; %Distance from centroid to (0,0,0)
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = datapoints(vertice,:) + dcm;
end
end

function xyTranslate(datapoints, containerAlpha, xTheta, yTheta)
%Input: datapoints - written as column vectors [x y z];
%Translates datapoints on xy-plane in theta direction until aggregates' max(x) or max(y) reaches max(x) or
%max(y)of container.
    step = 0
    while true
        xTranslate = 0.25*cos(xTheta);
        yTranslate = 0.25*sin(yTheta);
        for i = 1:length(datapoints(:,1))
                datapoints(i,1) = datapoints(i,1) + xTranslate;
                datapoints(i,2) = datapoints(i,2) + yTranslate;
        end
        pointCheck = inShape(containerAlpha, datapoints);
        if sum(pointCheck) > 0
            for i = 1:length(datapoints(:,1))
                datapoints(i,1) = datapoints(i,1) - xTranslate;
                datapoints(i,2) = datapoints(i,2) - yTranslate;
            end
            disp("Translation done")
            break
        end
        step = step + 1
    end
end

