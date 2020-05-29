clc
clear

%% Generating Repository of Aggregates
%Load in all of the aggregates
pullFiles
for i = 1:samples
    structname = strcat('agg',num2str(i));
    repos.(structname).Points = normalize(aggList{i}.Points);
    repos.(structname).ConnectivityList = aggList{i}.ConnectivityList;
    
    [Vbox, length, width, height] = getVolBox(aggList{i}.Points);
    repos.(structname).BoxHeight = height;
    repos.(structname).BoxWidth = width;
    repos.(structname).BoxLength = length;
    repos.(structname).Volume = getVolMesh(aggList{i}.Points);
    repos.(structname).VolumeFraction = getVolMesh(aggList{i}.Points)./Vbox;
end

[volumefraction, theta_x, theta_y, theta_z, datapoints] = optimize(repos.agg1.Points,repos.agg1.VolumeFraction);

%OUTPUT: A matrix with all of the aggregates ID, include properties: Mesh
%Volume, Box Volume (for packing), Width, Length, and Height (dim of box
%volume)

%% FUNCTION 1: Normalize Each Aggregate
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

%% FUNCTION 2: Generating the Volume of Each Aggregate
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
function Vmesh = getVolMesh(datapoints)
object = alphaShape(datapoints);
Vmesh = volume(object);
end

%% FUNCTION 3: Optimize Volume Fraction
function [volumefraction, theta_x, theta_y, theta_z, datapoints] = optimize(datapoints,volumefraction)
%starting condition
theta_x = 0;
theta_y = 0;
theta_z = 0;

angles = linspace(-pi/12,pi/12,50);
for a= 1:length(angles)
    for b = 1:length(angles)
        for g = 1:length(angles)
            %rotate aggregate
            rotx = [1 0 0 0; 0 cos(angles(a)) -sin(angles(a)) 0; 0 sin(angles(a)) cos(angles(a)) 0; 0 0 0 1];
            roty = [cos(angles(b)) 0 sin(angles(b)) 0; 0 1 0 0; -sin(angles(b)) 0 cos(angles(b)) 0; 0 0 0 1];
            rotz = [cos(angles(g)) -sin(angles(g)) 0 0; sin(angles(g)) cos(angles(g)) 0 0; 0 0 1 0; 0 0 0 1];
            R = rotx*roty*rotz;
            
            datapointsnew = datapoints;
            %obtain new volume fraction 
            for i = 1:length(datapoints(:,1))
                x = datapoints(i,1);
                y = datapoints(i,2);
                z = datapoints(i,3);
                
                p = [x; y; z; 1];
                pnew= R*p;
                datapointsnew(i,1) = pnew(1);
                datapointsnew(i,2) = pnew(2);
                datapointsnew(i,3) = pnew(3);
            end
            
            datapointsnew = normalize(datapointsnew);
            Vbox = getVolBox(datapointsnew);
            volumefraction_new = getVolMesh(datapointsnew)./Vbox;
            
            %save the new configuration if it is better than the last...
            if volumefraction_new > volumefraction
                volumefraction = volumefraction_new;
                theta_x = angles(a);
                theta_y = angles(b);
                theta_z = angles(g);
                datapoints = datapointsnew;
            end
        end
    end
end
end