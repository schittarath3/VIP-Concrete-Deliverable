clc
clear
%% Driving Code
tic

repos = generateRepo('STL Files\Aggregates Processed 3\*.stl', 27);         %Creating initial directory from stl files of aggregate

startCoords = [1 1 1];                                                     %Creating cubelets to store place aggregates
cubeSize = 654;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);              

load('grainsizeresults4.mat');
aggRepo = generateDistributedRepo(repos, 3000, results, sieveSz);
[insertRepo, aggRepo] = insertAgg(aggRepo, testCubes, 27);  %Inserting aggregates into cubes and generating a new repo

%using tangent planes to pack larger aggregates
fn = fieldnames(aggRepo);
packedRepo = struct;
packRepo.(fn{1}) = aggRepo.(fn{1});
load('grownSpheres2.mat');

%%
maxIter = 5;
allRepos = cell(maxIter,1);
parfor repo = 1:maxIter,3;
    packRepo = tangentPlane(insertRepo); 
    sphereCell = VoidSearch(packRepo, 654, 350,  0.20, grownSpheres); %finding empty spaces
    sphAggRepo = Sphere2Agg(aggRepo, sphereCell, results, sieveSz, repo, 999); % associating aggregates with spheres
    aRepo = mergeRepos(sphAggRepo, packRepo);
    allRepos{repo} = aRepo;
end
   
mergedRepo = struct;
for repo = 1:maxIter
    mergedRepo = mergeRepos(mergedRepo, allRepos{repo});
end
finalRepo = removeOverlaps(mergedRepo);
finalRepo = translateToPoint(finalRepo, 5, [654/2 654/2 654/2]);

%% 
maxIter = 3;
allRepos = cell(maxIter,1);
parfor repo = 1:maxIter,3
    sphereCell = VoidSearch(finalRepo, 654, 350,  0.20, grownSpheres); %finding empty spaces
    sphAggRepo = Sphere2Agg(aggRepo, sphereCell, results, sieveSz, repo, 999); % associating aggregates with spheres
    allRepos{repo} = sphAggRepo;
end

for repo = 1:maxIter
    finalRepo = mergeRepos(finalRepo, allRepos{repo});
end
finalRepo = removeOverlaps(finalRepo);
finalRepo = translateToPoint(finalRepo, 5, [654/2 654/2 654/2]);
finalRepo = removeOverlaps(finalRepo);
%% 

try
    mkdir 'STL Files' 'Aggregates Out'
end

repoToSTL(finalRepo, 2) %Create stl files in '.../STL Files/Aggregates Out'
disp('done')
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


