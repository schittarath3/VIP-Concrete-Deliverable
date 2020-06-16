startCoords = [1 1 1];
cubeSize = 150;
nDivisions = 3;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);
% newCubeSize = cubeSize - (cubeSize/nDivisions);
% newNDivisions = nDivisions - 1;
% a = (cubeSize/nDivisions)/2;
% newStartCoords = [a a a];
% moreCubes = genCublets(newCubeSize, newNDivisions, newStartCoords, 2);
% testCubes = vertcat(testCubes, moreCubes);

function cubeCoords = genCublets(cubeSize, nDivisions, startCoords, iteration)
%Generates a n*n*n size (n=cubeLength) cube with x (x=nDivisions) divisions
%to insert aggregates within each cublets of ((cubeSize/nDivisions)^3)
%size
%Inputs: length of cube (cubeSize), number of divisions (nDivisions),
%the iteration (iteration == 1 or == 2);
%Output: a cell, cubeCoords, containing the coordinates of the verticies
%       for each mini cube
step = cubeSize/nDivisions;
nDivisions = nDivisions + 1;

cube = cell(nDivisions, nDivisions, nDivisions);
if iteration == 1
    for x= 1:nDivisions
        for y = 1:nDivisions
            for z = 1:nDivisions
                cube{x,y,z} = [(step*(x-1)) (step*(y-1)) (step*(z-1))];
            end
        end
    end
end
if iteration == 2
    for x= 1:nDivisions
        for y = 1:nDivisions
            for z = 1:nDivisions
                cube{x,y,z} = [(step*(x-1)+startCoords(1)) (step*(y-1)+startCoords(2))...
                    (step*(z-1)+startCoords(3))];
            end
        end
    end
end

num = 1;
nDivisionsA = nDivisions-1;
preCell = cell(nDivisionsA^3,1);
for x = 1:nDivisionsA
    for y = 1:nDivisionsA
        for z = 1:nDivisionsA
            coordCell = [cube{x,y,z}; cube{x,y+1,z};...
                                  cube{x+1,y,z}; cube{x+1,y+1,z};...
                                  cube{x,y,z+1}; cube{x,y+1,z+1};...
                                  cube{x+1,y,z+1}; cube{x+1,y+1,z+1}];
             preCell{num,1} = coordCell;
             num = num + 1;
        end
    end
end

cubeCoords = preCell
end