clc
clear
%% Driving Code
tic

%Creating initial directory from stl files of aggregate
repos = generateRepo('STL Files\Aggregates Processed 3\*.stl', 27);         

%Creating cubelets to store place aggregates
startCoords = [1 1 1];                                                     
cubeSize = 654;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);              

%Inserting the 27 largest aggregates into cubes and generating a new repo
load('grainsizeresults4.mat');
aggRepo = generateDistributedRepo(repos, 3000, results, sieveSz);
[insertRepo, aggRepo] = insertAgg(aggRepo, testCubes, 27);  

%using tangent planes to pack larger aggregates
fn = fieldnames(aggRepo);
packedRepo = struct;
packRepo.(fn{1}) = aggRepo.(fn{1});
load('grownSpheres2.mat');

%% Aggregate Distribution Packing Run 1
maxIter = 5;
allRepos = cell(maxIter,1);

%We create 5 seperate parfor loops to greedily find empty spaces using 
%spheres
parfor repo = 1:maxIter,3;
    packRepo = tangentPlane(insertRepo); 
   
    %finding empty spaces
    sphereCell = VoidSearch(packRepo, 654, 350,  0.20, grownSpheres); 
    
    % associating aggregates with spheres based on grain-size distribution
    sphAggRepo = Sphere2Agg(aggRepo, sphereCell, results, sieveSz, repo, 999); 
    
    %merge original repo from tangentPlane and newly created repo
    aRepo = mergeRepos(sphAggRepo, packRepo); 
    allRepos{repo} = aRepo;
end
   
%merge all 5 repos
mergedRepo = struct;
for repo = 1:maxIter
    mergedRepo = mergeRepos(mergedRepo, allRepos{repo});
end

%remove any overlaps using alphaShapes
finalRepo = removeOverlaps(mergedRepo);

%translate alll cubes to center of the cube
finalRepo = translateToPoint(finalRepo, 5, [654/2 654/2 654/2]);

%% Aggregate Distribution Packing Run 2

%Create 5 seperate parfor loops to greedily find empty spaces using 
%spheres
maxIter = 3;
allRepos = cell(maxIter,1);
parfor repo = 1:maxIter,3
    
     %finding empty spaces
    sphereCell = VoidSearch(finalRepo, 654, 350,  0.20, grownSpheres);
    
     % associating aggregates with spheres
    sphAggRepo = Sphere2Agg(aggRepo, sphereCell, results, sieveSz, repo, 999);
    allRepos{repo} = sphAggRepo;
end

%merge all 5 repos and remove any overlaps
for repo = 1:maxIter
    finalRepo = mergeRepos(finalRepo, allRepos{repo});
end

%remove any overlaps using alphaShapes
finalRepo = removeOverlaps(finalRepo);

%translate alll cubes to center of the cube
finalRepo = translateToPoint(finalRepo, 5, [654/2 654/2 654/2]);

%Final cleanup
finalRepo = removeOverlaps(finalRepo);
%% Output aggregate STL

%Create directory in folder
try
    mkdir 'STL Files' 'Aggregates Out'
end

%Create stl files in '.../STL Files/Aggregates Out' as
%All.stl as the main file and individual aggregates
repoToSTL(finalRepo, 2) 
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


