clear
clc
%% Driving Code
repos = generateRepo('STL Files\Aggregates Processed 3\*.stl', 9);         %Creating initial directory from stl files of aggregate

startCoords = [1 1 1];                                                     %Creating cubelets to store place aggregates
cubeSize = 561;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);              

aggRepo = insertAgg(repos, testCubes, 1, 3);                               %Inserting aggregates into cubes and generating a new repo

%Tangent function call goes here

repoToSTL(aggRepo);                                                        %Converting points into Repo into STL then plotting STL
plotSTL('STL Files\Aggregates Out');