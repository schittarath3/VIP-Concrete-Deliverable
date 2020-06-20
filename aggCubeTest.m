clc 
clear

fileList = dir('STL Files\Aggregates Processed 3\*.stl');
for ag = 1:length(fileList)
    try
    filename = fileList(ag,1).name;
    folderName = fileList(ag,1).folder;
    fileLoc = strcat(folderName, '\', filename);
    agg = stlread(fileLoc);
    
    %generating the points and connectivity
    pts = normalize(agg.Points);
    cnt = agg.ConnectivityList;

    %finding volume of mesh
    set(0,'DefaultFigureVisible','off')
    model = createpde;
    importGeometry(model,fileLoc);
    mesh = generateMesh(model);
    Vmesh = volume(mesh);
    figure
    [V, nf, nv] = Volume(pts,cnt,Vmesh,fileLoc,false);
    repos.(filename(1:end-4)).Vertices = nv;
    repos.(filename(1:end-4)).Faces = nf;
    
    %rotating each of the aggregates for a set orientation...
    angles = linspace(-pi/8,pi/8,5);
    tz = 0;
    for ty = 1:length(angles)
        for tx = 1:length(angles)
            for tz = 1:length(angles)
                nv = Rotate(nv,angles(tx),angles(ty),angles(tz));
                orientation = strcat('tx_indx',num2str(tx),'ty_indx',num2str(ty),'tz_indx',num2str(tz));
                repos.(filename(1:end-4)).Orientation.(orientation) = nv;
            end
        end
    end
    
    catch
        disp(['error with' + filename]);
    end
end

startCoords = [1 1 1];
cubeSize = 561;
nDivisions = 4;
testCubes = genCublets(cubeSize, nDivisions, startCoords, 1);
newCubeSize = cubeSize - (cubeSize/nDivisions);
newNDivisions = nDivisions - 1;
a = (cubeSize/nDivisions)/2;
newStartCoords = [a a a];
moreCubes = genCublets(newCubeSize, newNDivisions, newStartCoords, 2);
testCubes = vertcat(testCubes, moreCubes);

aggRepo = insertAgg(repos, testCubes, 1, 7);
% myRepo = growAgg(aggRepo, testCubes, 1.005);

%% Functions
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

cubeCoords = preCell;
end

function nv = Rotate(pts,tx,ty,tz)
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end

function [aggVol, nf, nv] = Volume(pts,cnt,meshv,filename,plot)
%calculating percentage based on face
set(0,'DefaultFigureVisible','off') %turn off any figures

mesh = trimesh(cnt,pts(:,1),pts(:,2),pts(:,3));
tol = 200;
numFaces = length(mesh.Faces(:,1));
redper = tol/numFaces;
[nf, nv] = reducepatch(mesh,redper);

%calculating shrink percentage
stepsz = .02;
aggVol = abs(stlVolume(nv',nf'));
aggVf = meshv./aggVol;

volfract = .80; %ideal vol fraction
while aggVf > volfract
    nvnew = nv*(1.0+stepsz); %gradually increasing representative shape
    aggVolnew = abs(stlVolume(nvnew',nf'));
    aggVfnew = meshv./aggVolnew;
    
    if aggVfnew < volfract
        break
    else
        nv = nvnew;
        aggVf = aggVfnew;
        stepsz = stepsz*.87;
    end
end

if plot == true
    set(0,'DefaultFigureVisible','on')
    TR = trimesh(cnt,pts(:,1),pts(:,2),pts(:,3));
    hold on
    TR.FaceAlpha = 1;
    TR.EdgeColor = 'b';
    
    TRred = trimesh(nf,nv(:,1),nv(:,2),nv(:,3));
    TRred.FaceAlpha = .9;
    TRred.EdgeColor = 'k';
    xlabel('X axis')
    ylabel('Y axis')
    zlabel('Z axis')
    
    title([filename ': Cell with Volume Fraction: ' num2str(aggVf)]);
    axis equal
    hold off
end
end

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
