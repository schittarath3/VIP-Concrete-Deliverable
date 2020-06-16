clc
clear

function plotaggregate(agg2,agg1,rotangle,NormVec,option)
%for visualization purposes...
pts2 = normalize(stlread(agg2).Points);
if ~isempty(rotangle)
pts2 = Rotate(pts2,0,rotangle,0);
end
pts1 = normalize(stlread(agg1).Points);

if ~isempty(NormVec)
pts2 = Translate(pts2,pts1,NormVec,option);
end
cnt = stlread(agg2).ConnectivityList;
trimesh(cnt,pts2(:,1),pts2(:,2),pts2(:,3),'FaceAlpha',0,'EdgeColor','k');
hold on
end

function result = optAngle(agg1,agg2,restrictface1,restrictface2)
%find the optimal angle that will allow the two normal vectors to have a 0
%degree between them (meaning perfectly parallel)
%Inputs: agg1 - the stationary/fixed aggregate
%             agg2 - the aggregate that needs to be "attached"/shifted
%Outputs: result - the optimal angle that will generate two parallel planes
a = linspace(0,pi./2);
optimalangle = [];
plane1 = tangentPlane(agg1,restrictface1,0,false);

parfor ang = 1:length(a)
    plane2 = tangentPlane(agg2,restrictface2,a(ang),false);
    radt = anglebwPlanes(plane1,plane2); %find the angle
    optimalangle = [optimalangle; radt]; %store the angle
end
[~,idx] = min(optimalangle); %after running through list of possibilities,
%determine the angle that generates the smallest theta 
result = a(idx);
end

function result = anglebwPlanes(normvec1,normvec2)
%find the angle between two normal vectors from two plane equations...
%Inputs:    normvec 1 - the normal vector of the first plane (simply give
%                 1x4 matrix of the generated equation from f'n tangentPlane
%                normvec2 - the normal vector of the second plane
%Outputs: 
magNormVec1 = norm(normvec1(1:3));
magNormVec2 = norm(normvec2(1:3));
%applying dot product to find the angle between two vectors...
result = acos(dot(normvec1(1:3),normvec2(1:3))./(magNormVec1.*magNormVec2));
end

function result = tangentPlane(filename,option,angle,figurecheck)
%to generate the equation of the tangential plane to a selected region of
%the aggregate (top of the aggregate, bottom, left, right). a visualization
%of the plane is optional.
%Inputs:    pts - the points of aggregates in nx3 matrix (x-,y-,z-
%                coordinates)
%               option - the selected region to find the plane of the
%               aggregate
%               angle - the angle the aggregate rotate by the y-axis 
%               figurecheck - enter 'true' to receive a visualize figure of
%               the tangential plane (else enter 'false' for calculations)
%Outputs: result - the equation of the plane written as a 1x4 matrix in [A
%              B C D] (reference to equation Ax+By+Cz+D=0 for a plane)

if contains(option,'left')
    split1 = 1; %axis where the aggregate is cut into halves
    split2 = 3; %axis where the aggregate is cut into halves again
    split3 = 2;
    restrictface = 1;
elseif contains(option,'right')
    split1 = 1;
    split2 = 3;
    split3 = 2;
    restrictface = 4;
elseif contains(option,'bottom')
    split1 = 3;
    split2 = 2;
    split3 = 1;
    restrictface = 1;
elseif contains(option,'top')
    split1 = 3;
    split2 = 2;
    split3 = 1;
    restrictface = 4;
end

pts1 = normalize(stlread(filename).Points);
pts1 = Rotate(pts1,0,angle,0);

%generating the tangent line of the top...
[~,idx] = sort(pts1(:,split1));
sortedpts1 = pts1(idx,:);
negidx = find(sortedpts1(:,split1) <= 0); 
if restrictface ==1 || restrictface == 2 
    newpts1 = sortedpts1(1:negidx(end,1),:); 
else
    newpts1 = sortedpts1(negidx(end,1)+1:end,:);
end

%finding if there are any flat faces
lowz = unique(round(newpts1(:,split1),3));
lowzidx = find(newpts1(:,split1) > lowz(1));

%removing flat face to generate representative plane...
%rotating the aggregate until the min x is more than 0...
flatface = newpts1(1:lowzidx-1,:);

if length(flatface) > 100
dTheta = 1; 
stepSz = .25;
if ~isempty(flatface) %there is flat face
    while min(flatface(:,split2)) < 0
        dTheta = dTheta + stepSz;
        flatface = Rotate(flatface,0,0,-pi./dTheta);
    end
    pts1 = Rotate(pts1,0,0,-pi/dTheta);
    newpts1 = Rotate(newpts1,0,0,-pi/dTheta);
end
    flat = true;
else
    flat = false;
end

if figurecheck == true
    figure
    cnt1 = stlread(filename).ConnectivityList;
    mesh1 = trimesh(cnt1,pts1(:,1),pts1(:,2),pts1(:,3));
    mesh1.FaceAlpha = 1;
    mesh1.EdgeColor = 'b';
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    hold on
end

%which tangent plane needs to be find depends on the region...
if flat == true && restrictface == 2
    newpts1 = flatface;
else %further restrict the sliced section...
    [~,idx] = sort(newpts1(:,split2));
    sortedpts2 = newpts1(idx,:);
    negidx = find(sortedpts2(:,split2) <= 0);
        if restrictface == 1 || restrictface == 3
        newpts1 = sortedpts2(1:negidx(end,1)-1,:);
        else
        newpts1 = sortedpts2(negidx(end,1)+1:end,:);
        end
end
if figurecheck == true
    plot3(newpts1(:,1),newpts1(:,2),newpts1(:,3),'.r')
    hold on
end

%calculating the centroid
x = newpts1(:,1); y = newpts1(:,2); z = newpts1(:,3);
xcm = sum(x)./length(x); ycm = sum(y)./length(y); zcm = sum(z)./length(z);
r0 = [xcm ycm zcm];

%finding closest point to centroid...
order = abs(newpts1(:,1:2)-r0(1:2));
[~,minidx] = min(order(:,1) + order(:,2));
r0 = newpts1(minidx,:);
% plot3(r0(1),r0(2),r0(3),'.','MarkerSize',50), remove '%' if you need to
% visualize the centroid of the sliced section
% hold on

%writing the equation for tangent plane...the objective is to obtain a list
%of localized points about the proposed "centroid" of the sliced section
%and obtain the vectors that describes the span of the x- and
%y-coordinates
ptsX = [];
for i = 1:length(newpts1)
    difX= newpts1(i,split2) - r0(split2);
    if abs(difX) < 15
        term = newpts1(i,:);
        ptsX = [ptsX; term];
    end
end
%obtaining the vector with the most difference in Y
[~,maxY] = max(ptsX(:,split3));
[~,minY] = min(ptsX(:,split3));
WidV = ptsX(maxY,:) - ptsX(minY,:);

%finding the point in the same y-axis (relatively)
ptsY = [];
for i = 1:length(newpts1)
    difY= newpts1(i,1) - r0(split3);
    if abs(difY) < 15
        term = newpts1(i,:);
        ptsY = [ptsY; term];
    end
end
DepV = ptsY(end,:) - ptsY(split2,:);

Normvec = cross(DepV,WidV); %crossing the two vectors to get the normal vector...
Normvec = Normvec./norm(Normvec); 

%limit x and y to only the max width/length (for visual purposes with the
%normal plane)
limxy = abs([max(newpts1(:,1)), max(newpts1(:,2)), min(newpts1(:,1)), min(newpts1(:,2))]);
limxy = max(limxy);

[X,Y] = meshgrid(-limxy:2:limxy);
%move r0 so that the z is now zmin...
if restrictface == 1 || restrictface == 2
    [~,minZ] = min(newpts1(:,split1));
    r0 = newpts1(minZ,:); %readjust to include the minimal point in plane
else
    [~,maxZ] = max(newpts1(:,split1));
    r0 = newpts1(maxZ,:);
end
D = -1*(Normvec(1)*r0(1) + Normvec(2)*r0(2) + Normvec(3)*r0(3)); %D of the plane equation...
Z = -1/Normvec(3)*(Normvec(1)*X + Normvec(2)*Y + D); %equation for the plane described in Z(x,y)
if figurecheck == true
surf(X,Y,Z,'FaceColor','w');
hold off
axis equal
end

result = [Normvec, D];
end

function [datapointsn, centroid] = normalize(datapoints)
%obtain the new datapoints of the x-,y-,z- coordinates of the aggregates
%centered at the origin..
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

function datapointsn = Translate(datapoints1,datapoints2,Normvec,option)
if contains(option,'left') 
    scalar = -1;
elseif contains(option,'right')
    scalar = 1;
elseif contains(option,'bottom')
    scalar = -1;
elseif contains(option,'top')
    scalar = 1;
end

agg2 = alphaShape(datapoints2); %stationary aggregate (reference geometry)

alpha = 0;
StepSz = 1; %how many times to multiply the normal vector (scalar multiplier...)
datapointsn = datapoints1; %create a new matrix to store potential coordinates 

locate = sum(inShape(agg2,datapoints1)); %determine if the two shapes are overlapping or not...
while locate~=0 %if the aggregates are not overlapping...(overlapping -> sum = 0)
    Normvecn = scalar*Normvec.*(alpha + StepSz);

    datapointsn(:,1) = datapointsn(:,1) + Normvecn(1);
    datapointsn(:,2) = datapointsn(:,2) + Normvecn(2);
    datapointsn(:,3) = datapointsn(:,3) + Normvecn(3);
    locate = sum(inShape(agg2,datapointsn));

    if locate ~=0 
        datapointsn = datapoints1; %reset
        StepSz = StepSz + 1;
    else
        break
    end
end
end

function nv = Rotate(pts,tx,ty,tz)
%rotate the aggregate using 3D rotational matrices
%inputs: pts - a nx3 matrix with the columns as x-,y-,z- coordinates 
%             tx - the degrees (in RAD) of the rotation about the x-axis
%             ty - the degrees (in RAD) of the rotation about the y-axis
%             tz - the degrees (in RAD) of the rotation about the z-axis
%outputs: nv - a nx3 matrix with the columns as x-,y-,z- coordinates
%               of the new coordinates for the rotated aggregate
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;
nv = rotm*pts';
nv = nv';
end