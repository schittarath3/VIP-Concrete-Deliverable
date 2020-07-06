function repos = generateRepo(folder, numSamples)
fileList = dir(folder);
samples = size(fileList);
if numSamples <= samples
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
        repos.(filename(1:end-4)).Vertices = nv;
        repos.(filename(1:end-4)).Faces = nf;
        repos.(filename(1:end-4)).OriginalPoints = pts;
        repos.(filename(1:end-4)).OriginalFaces = agg.ConnectivityList;

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
else 
    disp("Number of samples less than number in directory")
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