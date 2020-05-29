clc
clear

%% Generating Repository of Aggregates
%Load in all of the aggregates and shrink them by percentage (set this as a
%variable so we can experiment later) - to shrink, use the scale matrix..
original = 20; %Number of aggregates that must be contained inside the box

%s = .90; %Shrink so it is .9 of original size

%Mesh Repository: [x y z connectivitylist volumefraction l w h];
numagg = 41; %Number of aggregates
mesh_opt = [];
for aggregate = 1:numagg
    %Define the volume here
    %Vmesh
    for alpha = linspace(-pi/12,pi/12)
        for beta = linspace(-pi/12,pi/12)
            for gamma = linspace(-pi/12,pi/12)
                %Normalize here...

                rotx = [1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
                roty = [cos(beta) 0 sin(beta) 0; 0 1 0 0; -sin(beta) 0 cos(beta) 0; 0 0 0 1];
                rotz = [cos(gamma) -sin(gamma) 0 0; sin(gamma) cos(gamma) 0 0; 0 0 1 0; 0 0 0 1];
                R = rotx*roty*rotz;

                s = 1; %Keep original shrunken size
                tx = 0; ty = 0; tz = 0; %Keep object at origin
                S = [s 0 0 0; 0 s 0 0; 0 0 s 0; 0 0 0 1];
                T = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];

                TR = R*S*T;

                %new matrix = old matrix (create new matrix)
                %for loop here, replace the coordinates with new coordinates
                %after transformation...
                %datapoints - the new set of coordinates
                [Vbox, l, w, h] = get.VolBox(datapoints);
                
                %The volume fraction (Volume of Mesh/Volume of Box)...
                volfrac = Vmesh/Vbox;
            end %stop for rotation angles in x
        end% stop for rotation angles in y
    end% stop for rotation angles in z
end% stop with all aggregates

%OUTPUT: A matrix with all of the aggregates ID, include properties: Mesh
%Volume, Box Volume (for packing), Width, Length, and Height (dim of box
%volume)

%% FUNCTION 1: Normalize Each Aggregate
function datapoints = get.Normalize(datapoints)
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

%% FUNCTION 2: Generating the Volume of Each Aggregate
function [Vbox, l, w, h] = get.VolBox(datapoints)
%Input: datapoints - written as column vectors [x y z];
%Output: V - the volume of the container/box
%            l, w, h - dimensions of the box (length, width, height)
w = abs(max(datapoints(:,1)) - min(datapoints(:,1)));
l = abs(max(datapoints(:,2)) - min(datapoints(:,2)));
h = abs(max(datapoints(:,3)) - min(datapoints(:,3));
Vbox = w.*l.*h;
end

function Vmesh = get.VolMesh(filename)
model = createpde;
importGeometry(model,filename);
mesh = generateMesh(model);
Vmesh = volume(mesh);
end

%% PACKING THE AGGREGATES