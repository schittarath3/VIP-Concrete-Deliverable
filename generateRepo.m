function repos = generateRepo(folder, numSamples)
%Generate a repository of the aggregates with the reduce mesh (to save
%calculation time) along with their connectivity list and the original
%points stored in a struct. The repository also includes the points of the
%aggregates rotated at a random angles.
%Inputs:
%   folder - folder location name (string)
%   numSamples - number of aggregates desired

%Read the stl files inside of the folder...
fileList = dir(folder);
samples = length(fileList);
repos = struct;

for ag = 1:numSamples
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

    %Generating the repository of the aggregates...
    repos.(filename(1:end-4)).Vertices = nv;
    repos.(filename(1:end-4)).Faces = nf;
    repos.(filename(1:end-4)).OriginalPoints = pts;
    repos.(filename(1:end-4)).OriginalFaces = agg.ConnectivityList;

    %Rotating each of the aggregates for a set orientation...
    angles = linspace(-pi/8,pi/8,5);
    tz = 0;
    for ty = 1:length(angles)
        for tx = 1:length(angles)
            for tz = 1:length(angles)
                nv = Rotate(nv,angles(tx),angles(ty),angles(tz));
                orientation = strcat('tx_indx',num2str(tx),'ty_indx',num2str(ty),'tz_indx',num2str(tz));
                repos.(filename(1:end-4)).Orientation.(orientation) = nv;

                %finding the maximum lengths
                repos.(filename(1:end-4)).xLength = max(nv(:,1)) - min(nv(:,1));
                repos.(filename(1:end-4)).yLength = max(nv(:,2)) - min(nv(:,2));
                repos.(filename(1:end-4)).zLength = max(nv(:,3)) - min(nv(:,3));
            end
        end
    end
    
    catch
        disp(['error with ' filename]);
    end
end
end
%% Functions
function nv = Rotate(pts,tx,ty,tz)
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end

function [aggVol, nf, nv] = Volume(pts,cnt,meshv,filename,plot)
%Finding the volume of the aggregates along with the connectivity list and
%the new coordinates of the reduced mes.

%Calculating percentage based on the number of faces...
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

volfract = .95; %Ideal volume fraction, this is the finite number that will stop
%the aggregate from further reducing. If the mesh is too reduced, it will
%be inside of the original and be too small.
while aggVf > volfract
    nvnew = nv*(1.0+stepsz); %Gradually increasing representative shape
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

if plot == true %If a plot is needed to see how the new reduced mesh compares with the original mesh.
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