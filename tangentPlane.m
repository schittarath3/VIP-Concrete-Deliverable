function aggRepov2 = tangentPlane(repos)
%Generate a repository with the coordinates of the packed aggregates by
%calculating the normal vectors of the tangent plane for a selected face of
%the aggregates.
%Inputs:
%   repos - the repository of the aggregates in appropriate cube cells
%   containing their position number, coordinates, and the connectivity
%   list generated from read stl. 
%Outputs:
 %  aggRepov2 - the new repository of the aggregates with the new
 %  coordinates of the packed aggregates.
fields = fieldnames(repos);

%Sorting the given repo in order of cube number...
newfields = fields;
for idx = 1:length(fields)
    orderidx = repos.(fields{idx}).cubeNum;
    newfields(orderidx,1) = fields(idx); 
end
fields = newfields;

%Defining the center aggregate where successive aggregates will be added by
%translating towards it (this is defined as the 14th position in the cube
%cells).
figure(1)
ang = linspace(-pi/8,pi/8,5);
rotOrient = repos.(fields{14}).Orientation;
pts14 = Rotate(repos.(fields{14}).OriginalPoints,ang(rotOrient(1)),ang(rotOrient(2)),ang(rotOrient(3)));
[~,cm14] = normalize(repos.(fields{14}).Points);
pts14 = pts14 + cm14; aggpts = pts14;
cnt14 = repos.(fields{14}).OriginalFaces;
trimesh(cnt14,pts14(:,1),pts14(:,2),pts14(:,3),'EdgeColor','blue');

aggRepov2.(fields{14}).Faces = cnt14;
aggRepov2.(fields{14}).Points = pts14;
hold on

sequence = randperm(27);
for idx = sequence
    if idx~=14 
        %Rotating the aggregates according the correct orientation 
        rotaOrient = repos.(fields{idx}).Orientation;
        pts = Rotate(repos.(fields{idx}).OriginalPoints,ang(rotaOrient(1)),ang(rotaOrient(2)),ang(rotaOrient(3)));
    
        %Translate the added aggregates towards the center aggregate,
        %making sure that the respective tangent planes are aligned at
        %nearly 0 degrees angle.
        [~,cm] = normalize(repos.(fields{idx}).Points); pts = pts+cm;
        cnt = repos.(fields{idx}).OriginalFaces;
        
        pts = translate2Pt(repos.(fields{idx}).cubeNum,pts,aggpts,pts14,clustercm,cm);
        aggpts = [aggpts; pts];

        trimesh(cnt,pts(:,1),pts(:,2),pts(:,3),'EdgeColor','blue');
        aggRepov2.(fields{idx}).Faces = cnt;
        aggRepov2.(fields{idx}).Points = pts;
        hold on
        xlabel('x'); ylabel('y'); zlabel('z');
        axis equal;
    end
end

%Defining the maximum dimension of the imaginary box that the aggregates
%will be packed into. Any aggregates that are translated and outside of the
%box will be discarded.*
marx = 10; mary = 5; marz = -10;
maxZ = max(aggpts(:,3))+marz; minZ = min(aggpts(:,3))-marz;
maxX = max(aggpts(:,1))+marx; minX = min(aggpts(:,1))-marx;
maxY = max(aggpts(:,2))+mary; minY = min(aggpts(:,2))-mary;

%Algorithm for additional aggregates (past the original 27 assigned inside
%of the generated cube).
for addAgg = 28:length(fields)
    rotaOrient = repos.(fields{addAgg}).Orientation;
    pts = normalize(Rotate(repos.(fields{addAgg}).OriginalPoints,ang(rotaOrient(1)),ang(rotaOrient(2)),ang(rotaOrient(3))));
    [~,cm] = normalize(repos.(fields{addAgg}).Points); 
    
    %Randomly rotate the aggregate by a small angle
    pts = pts+Rotate(cm,pi*(1/ randi([5 15],1)),pi*(1/ randi([5 15],1)),pi*(1/randi([5 15],1)));
    pts = translate2Pt(repos.(fields{addAgg}).cubeNum,pts,aggpts,cm14,cm);

    %Checking if the aggregates are translated inside of the box. Any
    %aggregates that are too big for the defined box will be removed.
    if max(pts(:,1)) > maxX || min(pts(:,1)) < minX || max(pts(:,2)) > maxY || min(pts(:,2)) < minY || max(pts(:,3)) > maxZ || max(pts(:,3)) < minZ
        continue
    else
        aggpts = [aggpts; pts];
        cnt = repos.(fields{addAgg}).OriginalFaces;
        aggRepov2.(fields{addAgg}).Faces = cnt;
        aggRepov2.(fields{addAgg}).Points = pts;
        trimesh(cnt,pts(:,1),pts(:,2),pts(:,3),'EdgeColor','green');
        hold on
    end
end
hold off
end

%% FUNCTION CODE
function agg = translate2Pt(cubeidx,agg,cluster,centerpts,clustercm,aggcm)
%Given two aggregates, the aggregates will be translated with one fixed and
%another towards its center. The translated aggregate will attempt to
%rotate so that the tangent planes of the two faces are 0-degrees. 
%Inputs:
%   cubeidx - the cube position of the aggregate
%   agg - the translated aggregate coordinates with the columns
%   representing the (x,y,z) coordinates
%   cluster - the fixed aggregate or cluster of aggregate with the columns
%   representing the (x,y,z) coordinates
%   clustercm - the center point of the aggregate/cluster
%   aggcm - the center point of the translated aggregate

%The selected face of the aggregates depending on their position. This face
%tells the code which face it will calculate the normal vector of.
switch cubeidx
    case {2, 11, 20}
        restrictface1 = [2 3 1 1]; rotidx = 2;
        if cubeidx == 2
            restrictface2 = [1 2 3 4];
        elseif cubeidx == 11
            restrictface2 = [2 3 1 3];
        elseif cubeidx == 20
            restrictface2 = [1 2 3 2];
        end
        
    case {8, 17, 26}
        restrictface1 = [2 1 3 3]; rotidx = 2;
        if cubeidx == 8
            restrictface2 = [1 2 3 3];
        elseif cubeidx == 17
            restrictface2 = [2 3 1 1];
        elseif cubeidx == 26
            restrictface2 = [1 3 2 1];
        end
        
    case 13
        restrictface1 = [3 2 1 1];
        restrictface2 = [3 2 1 4];
        rotidx = 3;
        
    case 15
        restrictface1 = [3 2 1 4];
        restrictface2 = [3 2 1 1];
        rotidx = 3;
        
    case {1, 4, 7, 10, 16, 19, 22, 25}
        restrictface1 = [3 1 2 1]; rotidx = 2;
        if cubeidx == 1 || cubeidx == 4 || cubeidx == 7
            restrictface2 = [1 3 2 3];
        elseif cubeidx == 19 || cubeidx == 22 || cubeidx == 25
            restrictface2 = [1 3 2 2];
        elseif cubeidx == 16
            restrictface2 = [2 1 3 2];
        elseif cubeidx == 10
            restrictface2 = [2 1 3 3];
        end
        
    case {3, 6, 9, 12, 18, 21, 24, 27}
        restrictface1 = [3 2 1 3]; rotidx = 2;
        if cubeidx == 3 || cubeidx == 6 || cubeidx == 9
            restrictface2 = [1 3 2 4]';
        elseif cubeidx == 12 || cubeidx == 18 || cubeidx == 21
            restrictface2 = [2 1 3 3];
        elseif cubeidx == 24
            restrictface2 = [2 3 1 1];
        elseif cubeidx == 27
            restrictface2 = [2 3 1 2];
        end
        
    case 23 
        restrictface1 = [1 3 2 3];
        restrictface2 = [1 2 3 1];
        rotidx = 2;
    case 5 
        restrictface1 = [1 2 3 1];
        restrictface2 = [1 3 2 3];
        rotidx = 2;
end

%Finding if the translated aggregate touches any other aggregate or has
%translated inside of the center aggregate using inShape.
clusterA = alphaShape(cluster);
locate = sum(inShape(clusterA,agg)); 

try
    result = optAngle(normalize(agg),centerpts,restrictface1,restrictface2,rotidx);
catch
    result = 0;
end
    if rotidx == 2
        agg = Rotate(normalize(agg),0,result,0) + aggcm;
    else
        agg = Rotate(normalize(agg),0,0,result) + aggcm;
    end

    while locate == 0 %0 means that no points are touching (none of the points are inside)
        [~,aggcm] = normalize(agg);
        tVect = (clustercm - aggcm)./75; %The translation vector defined as the direction between 
        %center points of the fixed aggregates/cluster and translating aggregate.
        agg = agg + tVect; %Slowly approaching towards the center
        locate = sum(inShape(clusterA,agg));
            if locate~=0
                agg = agg - tVect; %If it is inside, move back a step size to remain outside.
            end
    end
end

 function result = optAngle(pts1,pts2,restrictface1,restrictface2,opt)
%Find the optimal angle to rotate the translating aggregate so that the two
%calculated tangent planes are at 0-degrees with each other.
a = linspace(-pi/3,pi/3);
optimalangle = [];
plane1 = tangentP(pts1,restrictface1,[0 0 0],false);

parfor ang = 1:length(a) %Testing each rotation to find the optimal.
    if opt == 2 %y-axis
        rots = [0, a(ang), 0];
    else 
        rots = [0, 0, a(ang)];
    end
    plane2 = tangentP(pts2,restrictface2,rots,false);
    radt = anglebwPlanes(plane1,plane2);
    optimalangle = [optimalangle; radt];
end

[~,idx] = min(optimalangle);
result = a(idx);
end

function result = anglebwPlanes(normvec1,normvec2)
%Calculate the angle between two tangent plane using dot product.
magNormVec1 = norm(normvec1(1:3));
magNormVec2 = norm(normvec2(1:3));

result = acos(dot(normvec1(1:3),normvec2(1:3))./(magNormVec1.*magNormVec2));
end

function result = tangentP(pts1,option,angle,figurecheck)
%Calculate the normal vector of the tangent plane of the selected face of
%the aggregate. To isolate the face, the object is split by the three axes
%(x,y,z) and the face is defined by the negative or positive of the axes.
%Using those points, the extreme point is used as an "anchor" to generate
%the vectors that approximates the plane of the surface. 

split1 = option(1);
split2 = option(2);
split3 = option(3);
restrictface = option(4);
pts1 = Rotate(pts1,angle(1),angle(2),angle(3));

%Generating the tangent line of the top...
[~,idx] = sort(pts1(:,split1));
sortedpts1 = pts1(idx,:);
negidx = find(sortedpts1(:,split1) <= 0); 
if restrictface ==1 || restrictface == 2 
    newpts1 = sortedpts1(1:negidx(end,1),:); 
else
    newpts1 = sortedpts1(negidx(end,1)+1:end,:);
end

[~,idx] = sort(newpts1(:,split2));
sortedpts2 = newpts1(idx,:);
negidx = find(sortedpts2(:,split2) <= 0);
    if restrictface == 1 || restrictface == 3
    newpts1 = sortedpts2(1:negidx(end,1)-1,:);
    else
    newpts1 = sortedpts2(negidx(end,1)+1:end,:);
    end

if figurecheck == true
    plot3(newpts1(:,1),newpts1(:,2),newpts1(:,3),'.r')
    hold on
end

%calculating the centroid
x = newpts1(:,1);
y = newpts1(:,2);
z = newpts1(:,3);

xcm = sum(x)./length(x);
ycm = sum(y)./length(y);
zcm = sum(z)./length(z);
r0 = [xcm ycm zcm];

%Finding closest point to centroid...
order = abs(newpts1(:,1:2)-r0(1:2));
[~,minidx] = min(order(:,1) + order(:,2));
r0 = newpts1(minidx,:);
% plot3(r0(1),r0(2),r0(3),'.','MarkerSize',50)
% hold on

%Writing the equation for tangent plane...
%Finding the point in the same x-axis (relatively)
ptsX = [];
for i = 1:length(newpts1)
    difX= newpts1(i,split2) - r0(split2);
    if abs(difX) < 25
        term = newpts1(i,:);
        ptsX = [ptsX; term];
    end
end
%Obtaining the vector with the most difference in Y
[~,maxY] = max(ptsX(:,split3));
[~,minY] = min(ptsX(:,split3));
WidV = ptsX(maxY,:) - ptsX(minY,:);

%Finding the point in the same y-axis (relatively)
ptsY = [];
for i = 1:length(newpts1)
    difY= newpts1(i,1) - r0(split3);
    if abs(difY) < 25
        term = newpts1(i,:);
        ptsY = [ptsY; term];
    end
end
DepV = ptsY(end,:) - ptsY(split2,:);

Normvec = cross(DepV,WidV);
Normvec = Normvec./norm(Normvec);

%Limit x and y to only the max width/length
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
D = -1*(Normvec(1)*r0(1) + Normvec(2)*r0(2) + Normvec(3)*r0(3));
Z = -1/Normvec(3)*(Normvec(1)*X + Normvec(2)*Y + D);
if figurecheck == true
surf(X,Y,Z,'FaceColor','w');
hold off
axis equal
end

result = [Normvec, D];
end

function [datapointsn, centroid] = normalize(datapoints)
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
%Rotational matrix
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end