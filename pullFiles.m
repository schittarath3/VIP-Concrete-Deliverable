%comment
fileList = dir('STL Files\Aggregates Processed 3\*.stl')
fileList = struct2cell(fileList);
aggList = cell(1, length(fileList));
samples = size(fileList);
meshvol = zeros(1,samples(2));
for i = 1:samples(2)
    %obtain stl
    folder = fileList{2, i};
    file = fileList{1, i};
    aggName = file(1:end-4);
    fn = strcat(folder, '\', file);
    aggList{i} = stlread(fn);
   
    %generating repository of aggregates
    repos.(aggName).Points = normalize(aggList{i}.Points);
    repos.(aggName).ConnectivityList = aggList{i}.ConnectivityList;
    
    %finding volume fraction of boxes to pack
    model = createpde;
    importGeometry(model,fn);
    mesh = generateMesh(model);
    Vmesh = volume(mesh);

    [Vbox, length, width, height] = getVolBox(aggList{i}.Points);
    repos.(aggName).BoxHeight = height;
    repos.(aggName).BoxWidth = width;
    repos.(aggName).BoxLength = length;
    repos.(aggName).Volume = Vmesh;
    repos.(aggName).VolumeFraction = Vmesh/Vbox;
end

function [Vbox, l, w, h] = getVolBox(datapoints)
%Input: datapoints - written as column vectors [x y z];
%Output: V - the volume of the container/box
%            l, w, h - dimensions of the box (length, width, height)
datapoints = normalize(datapoints);

w = abs(max(datapoints(:,1)) - min(datapoints(:,1)));
l = abs(max(datapoints(:,2)) - min(datapoints(:,2)));
h = abs(max(datapoints(:,3)) - min(datapoints(:,3)));
Vbox = w.*l.*h;
end

function datapoints = normalize(datapoints)
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
dcm = centroid - [0 0 0]; %Distance from centroid to (0,0,0)
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = [datapoints(vertice,:) - dcm];
end
end