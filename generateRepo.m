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
repos = struct;

for ag = 1:numSamples
    try
    filename = fileList(ag,1).name;
    folderName = fileList(ag,1).folder;
    fileLoc = strcat(folderName, '\', filename);
    agg = stlread(fileLoc);

    %Generating the points and connectivity
    pts = normalize(agg.Points);
    cnt = agg.ConnectivityList;

    %Finding volume of mesh
    set(0,'DefaultFigureVisible','off')
    model = createpde;
    importGeometry(model,fileLoc);
    mesh = generateMesh(model);
    Vmesh = volume(mesh);
    figure
    [new_cnt, new_pts] = Volume(pts,cnt,Vmesh,fileLoc,false);

    %Generating the repository of the aggregates...
    repos.(filename(1:end-4)).Vertices = new_pts;
    repos.(filename(1:end-4)).Faces = new_cnt;
    repos.(filename(1:end-4)).OriginalPoints = pts;
    repos.(filename(1:end-4)).OriginalFaces = agg.ConnectivityList;

    %Rotating each of the aggregates for a set orientation...
    angles = linspace(-pi/8,pi/8,5);
    for ty = 1:length(angles)
        for tx = 1:length(angles)
            for tz = 1:length(angles)
                rotate_pts = Rotate(new_pts,angles(tx),angles(ty),angles(tz));
                orientation = strcat('tx_indx',num2str(tx),'ty_indx',num2str(ty),'tz_indx',num2str(tz));
                repos.(filename(1:end-4)).Orientation.(orientation) = rotate_pts;

                %Finding the maximum lengths
                repos.(filename(1:end-4)).xLength = max(new_pts(:,1)) - min(new_pts(:,1));
                repos.(filename(1:end-4)).yLength = max(new_pts(:,2)) - min(new_pts(:,2));
                repos.(filename(1:end-4)).zLength = max(new_pts(:,3)) - min(new_pts(:,3));
            end
        end
    end
    
    catch
        disp(['error with ' filename]);
    end
end
end

%% Functions

function [new_cnt, new_pts] = Volume(pts,cnt,mesh_volume,filename,plot)
%Finding the volume of the aggregates along with the connectivity list and
%the new coordinates of the reduced mesh.
%Inputs:
    %pts - coordinates of the aggregates with each column representing x,
    %y, z
    %cnt - connectivity list of the aggregate
    %mesh_volume - the volume of the mesh
    %filename - the string of the stl file name in folder
    %plot - Plot a figure of the reduced mesh conpared to the original
    %aggregates (true to plot, false to disable figure)
%Outputs:
    %new_cnt - new connectivity list of reduced mesh
    %new_pts - new aggregate coordinates of reduced mesh

%Calculating percentage based on the number of faces...
set(0,'DefaultFigureVisible','off') %turn off any figures
mesh = trimesh(cnt,pts(:,1),pts(:,2),pts(:,3));
tol = 250;
numFaces = length(mesh.Faces(:,1));
redper = tol/numFaces; 
[new_cnt, new_pts] = reducepatch(mesh,redper);

%Calculating shrink percentage
stepsz = .02;
aggVol = abs(stlVolume(new_pts',new_cnt'));
aggVf = mesh_volume./aggVol;

volfract = 1; %Ideal volume fraction, this is the finite number that will stop
%the aggregate from further reducing. If the mesh is too reduced, it will
%be inside of the original and be too small.
while aggVf > volfract
    iter = 0;
    nvnew = new_pts*(1.0+stepsz); %Gradually increasing representative shape
    aggVolnew = abs(stlVolume(nvnew',new_cnt'));
    aggVfnew = mesh_volume./aggVolnew;
    
    if aggVfnew < volfract
        break
    else
        iter = iter +1;
        new_pts = nvnew;
        aggVf = aggVfnew;
        stepsz = stepsz.*(1/iter);
    end
end

%If a plot is needed to see how the new reduced mesh compares with the original mesh.
if plot == true 
    set(0,'DefaultFigureVisible','on')
    TR = trimesh(cnt,pts(:,1),pts(:,2),pts(:,3));
    hold on
    TR.FaceAlpha = 1;
    TR.EdgeColor = 'b';
    
    TRred = trimesh(new_cnt,new_pts(:,1),new_pts(:,2),new_pts(:,3));
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
%Translate the aggregates to the origin (0,0,0)
%Inputs: 
    %datapoints - coordinates of the aggregates with each column
    %representing x,y,z
%Outputs:   
    %datapointsn - translated coordinates
    
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

function nv = Rotate(pts,tx,ty,tz)
%Rotate the aggregates using a rotational matrix
%Inputs:
    %pts - coordinates of the aggregates with each column representing
    %x,y,z
    %tx - angle about the x-axis
    %ty - angle about the y-axis
    %tz - angle about the z-axis
%Outputs:
    %nv - rotated coordinates
    
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end