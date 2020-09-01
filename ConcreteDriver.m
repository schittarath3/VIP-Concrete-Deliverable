clc
clear
%% Driving Code
tic

repos = generateRepo('STL Files\Aggregates Processed 3\*.stl', 27);         %Creating initial directory from stl files of aggregate

startCoords = [1 1 1];                                                     %Creating cubelets to store place aggregates
cubeSize = 654;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);              

load('grainsizeresults.mat');
aggRepo = generateDistributedRepo(repos, 1000, results, sieveSz);
[insertRepo, aggRepo] = insertAgg(aggRepo, testCubes, 27);  %Inserting aggregates into cubes and generating a new repo

firstPackRepo = tangentPlane(insertRepo); %using tangent planes to pack larger aggregates
% 
% spherePoints = stlread('sphere.stl').Points;
% sphereCell = VoidSearch(firstPackRepo, 654, 300,  0.70, spherePoints); %finding empty spaces
% [sphAggRepo, a] = Sphere2Agg(aggRepo, sphereCell, results, sieveSz, 999); % associating aggregates with spheres
% 
% mergedRepo = mergeRepos(sphAggRepo, firstPackRepo);
% 
% finalRepo = removeOverlaps(mergedRepo);
% 
% try
%     mkdir 'STL Files' 'Aggregates Out'
% end
% 
% repoToSTL(finalRepo, 2) %Create stl files in '.../STL Files/Aggregates Out'
% disp('done')
% toc
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


