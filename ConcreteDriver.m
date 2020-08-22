clc
clear
%% Driving Code
tic
repos = generateRepo('STL Files\Aggregates Processed 3\*.stl', 27);         %Creating initial directory from stl files of aggregate

startCoords = [1 1 1];                                                     %Creating cubelets to store place aggregates
cubeSize = 654;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);              

load('grainsizeresults.mat')
aggRepo = generateDistributedRepo(repos, 1000, results, sieveSz)
[insertRepo, aggRepo] = insertAgg(aggRepo, testCubes, 27);  %Inserting aggregates into cubes and generating a new repo

%Tangent function call goes here
finalRepo = tangentPlane(insertRepo);
[finalRate, totalVolume] = coverageRate(finalRepo, 218^3);
%finalRepo = shrinkByOrigin(finalRepo, .10769);

repoToSTL(finalRepo);                                                        %Converting points into Repo into STL then plotting STL
plotSTL('STL Files\Aggregates Out');
toc
%% Functions
function mergedRepo = mergeRepos(repoA,repoB)
%Merges two aggRepos
    mergedRepo = repoA;
    fn = fieldnames(repoB);
    structBLen = length(fn);
    for i = 1:structBLen
        mergedRepo.(fn{i}) = repoB.(fn{i});
    end
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
