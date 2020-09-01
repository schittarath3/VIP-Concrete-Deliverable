function aggRepo = tangentPlane(repos)
%Generate a repository with the coordinates of the packed aggregates by
%calculating the normal vectors of the tangent plane for a selected face of
%the aggregates.
%Inputs:
%   repos - the repository of the aggregates in appropriate cube cells
%   containing their position number, coordinates, and the connectivity
%   list generated from read stl. 
%Outputs:
 %  aggRepo - the new repository of the aggregates with the new
 %  coordinates of the packed aggregates.
fields = fieldnames(repos);

%Sorting the given repo in order of cube number...
newfields = fields;
for idx = 1:length(fields)
    orderidx = repos.(fields{idx}).cubeNum;
    newfields(orderidx,1) = fields(idx); 
end
fields = newfields;

%Zero matrix to store each of the aggregates center...
cm = zeros(27,3);
%Selection of angles
ang = linspace(-pi/8,pi/8,5);

%Defining the center aggregate where successive aggregates will be added by
%translating towards it (this is defined as the 14th position in the cube
%cells).
figure(1)
pts14 = repos.(fields{14}).OriginalPoints;
cnt14 = repos.(fields{14}).OriginalFaces;
[~,cm(14,1:3)] = normalize(repos.(fields{14}).Points);
pts14 = pts14 + cm(14); 
pts_index = ones(length(pts14),1)*14;
aggpts = [pts14, pts_index];

aggRepo.(fields{14}).OriginalFaces = cnt14;
aggRepo.(fields{14}).OriginalPoints = pts14;
aggRepo.(fields{14}).Points = repos.(fields{14}).Points;
aggRepo.(fields{14}).Faces = repos.(fields{14}).Faces;
aggRepo.(fields{14}).Orientation = repos.(fields{14}).Orientation;
aggRepo.(fields{14}).Diameter = repos.(fields{14}).Diameter;
aggRepo.(fields{14}).bin= repos.(fields{14}).bin;
aggRepo.(fields{14}).cubeNum= repos.(fields{14}).cubeNum;

trimesh(cnt14,pts14(:,1),pts14(:,2),pts14(:,3),'EdgeColor','blue');
hold on

%Sequence to pack each of the aggregates according to cube cell number...
sequence = [5,17,13,11,15,23,6,12,18,24,4,10,16,22,1,2,3,7,8,9,19,20,21,25,26,27];
for idx = 1:length(sequence)
    %Rotating the aggregates according the correct orientation 
    aggNum = sequence(idx);
    pts = repos.(fields{aggNum}).OriginalPoints;
    
    %Translate the added aggregates towards the center aggregate,
    %making sure that the respective tangent planes are aligned at
    %nearly 0 degrees angle.
    [~,cm(aggNum,1:3)] = normalize(repos.(fields{aggNum}).Points); pts = pts+cm(aggNum,1:3);
        
    if idx<=6
        [pts, result] = translate2Pt(repos.(fields{aggNum}).cubeNum,pts,aggpts,cm(aggNum,1:3),14,cm(14,1:3));
    elseif idx>6 && idx<=14
        [pts, result] = translate2Pt(repos.(fields{aggNum}).cubeNum,pts,aggpts,cm(aggNum,1:3),15,cm(14,1:3));
    elseif idx>14 && idx<=16
        [pts, result] = translate2Pt(repos.(fields{aggNum}).cubeNum,pts,aggpts,cm(aggNum,1:3),13,cm(14,1:3));
    else
        [pts, result] = translate2Pt(repos.(fields{aggNum}).cubeNum,pts,aggpts,cm(aggNum,1:3),[],cm(14,1:3));
    end
    cnt = repos.(fields{aggNum}).OriginalFaces;
    pts_index = ones(length(pts),1)*aggNum;
    newpts = [pts, pts_index];
    aggpts = [aggpts; newpts];

        aggRepo.(fields{aggNum}).OriginalPoints = pts;
        aggRepo.(fields{aggNum}).OriginalFaces = cnt;
        
        [~,finalcm] = normalize(pts);
        reduce_pts = Rotate(normalize(repos.(fields{aggNum}).Points),0,result(1),result(2)) + finalcm;
        aggRepo.(fields{aggNum}).Points = reduce_pts;
        aggRepo.(fields{aggNum}).Faces = repos.(fields{aggNum}).Faces;
        aggRepo.(fields{aggNum}).Orientation = repos.(fields{aggNum}).Orientation;
        aggRepo.(fields{aggNum}).Diameter = repos.(fields{aggNum}).Diameter;
        aggRepo.(fields{aggNum}).bin= repos.(fields{aggNum}).bin;
        aggRepo.(fields{aggNum}).cubeNum= repos.(fields{aggNum}).cubeNum;

        trimesh(cnt,newpts(:,1),newpts(:,2),newpts(:,3),'EdgeColor','red');
        hold on
        trimesh(repos.(fields{aggNum}).Faces,reduce_pts(:,1),reduce_pts(:,2),reduce_pts(:,3),'EdgeColor','blue');
        hold on
        xlabel('x'); ylabel('y'); zlabel('z');
        axis equal;
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
    pts = translate2Pt(repos.(fields{addAgg}).cubeNum,pts,aggpts,cm(aggNum,1:3),[],cm(14,1:3));

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

%% FUNCTIONS
function [agg_pts, agg_rotate] = translate2Pt(cubeidx,agg_pts,cluster_pts,agg_cm,adj_aggnum,target)
%Given two aggregates, the aggregates will be translated with one fixed and
%another towards its center. The translated aggregate will attempt to
%rotate so that the tangent planes of the two faces are 0-degrees. 
%Inputs:
%   cubeidx - the cube position of the aggregate
%   agg_pts - the translated aggregate coordinates with the columns
%   representing the (x,y,z) coordinates
%   cluster_pts - the fixed aggregate or cluster of aggregate with the columns
%   representing the (x,y,z) coordinates
%   agg_cm - the center of the aggregate
%   adj_aggnum - the aggregate's adjacent aggregate cube cell number/index 
%   target - the point the aggregate translates to
%Outputs:
%   agg_pts - the new translated aggregate coordinates with columns
%   representing the (x,y,z) coordinates
%   agg_rotate - the angle the aggregates rotated before translation/after
%   tangent plane calculations.

%The selected face of the aggregates depending on their position. This face
%tells the code which face it will calculate the normal vector of.
switch cubeidx
    case {5,11,13,15,17,23,6,12,18,24,4,10,16,22}
        if cubeidx == 4 || cubeidx == 5 || cubeidx == 6
            restrictface1 = [1 2 3 1];
            restrictface2 = [1 3 2 3];
        elseif cubeidx == 10 || cubeidx == 11 || cubeidx == 12
            restrictface1 = [1 2 3 1];
            restrictface2 = [1 3 2 3];
        elseif cubeidx == 13
            restrictface1 = [3 2 1 1];
            restrictface2 = [3 2 1 4];
        elseif cubeidx ==15
            restrictface1 = [3 2 1 4];
            restrictface2 = [3 2 1 1];
        elseif cubeidx == 16 || cubeidx == 17 || cubeidx == 18
            restrictface1 = [1 2 3 1];
            restrictface2 = [1 3 2 3];
        elseif cubeidx == 22 || cubeidx == 23 || cubeidx == 24
            restrictface1 = [1 3 2 3];
            restrictface2 = [1 2 3 1];
        end
        
        adj_agg = cluster_pts(find(cluster_pts(:,4) == adj_aggnum),1:3);
        agg_rotate = optAngle(agg_pts,adj_agg,restrictface1,restrictface2,[]);
        agg_pts = Rotate(normalize(agg_pts),0,agg_rotate(1),agg_rotate(2)) + agg_cm;
        
    case {3,9,21,27,1,7,19,25,2,8,20,26,27}
        %Calculate the angles between the faces of the adjacent aggregates
        %and rotate these aggregates appropriately so that there are spaces
        %between to rotate towards
            if cubeidx>=1 && cubeidx <=3
                adj_aggnum1 = cubeidx+3; 
                adj_aggnum2 = cubeidx+9;
                adj_restrictface1 = [2 3 1 4]; %front side
                adj_restrictface2 = [1 3 2 1]; %left side
                restrictface1 = [2 3 1 1]; %back side
                restrictface2 = [1 3 2 4]; %right side
            elseif cubeidx>=7 && cubeidx <=9
                adj_aggnum1 = cubeidx-3;
                adj_aggnum2 = cubeidx+9;
                adj_restrictface1 = [2 3 1 1]; %back side
                adj_restrictface2 = [1 3 2 1]; %left side
                restrictface1 = [2 3 1 4]; %front side
                restrictface2 = [1 3 2 4]; %right side
            elseif cubeidx>=19 && cubeidx <=21
                adj_aggnum1 = cubeidx-9;
                adj_aggnum2 = cubeidx+3;
                adj_restrictface1 = [1 3 2 4]; %right side
                adj_restrictface2 = [2 3 1 4]; %front side
                restrictface1 = [1 3 2 1]; %left side
                restrictface2 = [2 3 1 1]; %back side
            elseif cubeidx>=25 && cubeidx <=27
                adj_aggnum1 = cubeidx-3;
                adj_aggnum2 = cubeidx-9;
                adj_restrictface1 = [2 3 1 1]; %back side
                adj_restrictface2 = [1 3 2 4]; %right side
                restrictface1 = [2 3 1 4]; %front side
                restrictface2 = [1 3 2 1]; %left side
            end
        adj_agg1 = cluster_pts(find(cluster_pts(:,4) == adj_aggnum1),1:3);
        adj_agg2 = cluster_pts(find(cluster_pts(:,4) == adj_aggnum2),1:3);
        
        %Calculate the angle between the adjacent aggregates
        adj_ang1 = tangentP(adj_agg1,adj_restrictface1,[0 0 0],false);
        adj_ang2 = tangentP(adj_agg2,adj_restrictface2,[0 0 0],false);
        adj_ang = anglebwPlanes(adj_ang1,adj_ang2);
        
        %Calculate the angle between the normal vectors of the aggregate of
        %interest:
        agg_rotate = optAngle(agg_pts,agg_pts,restrictface1,restrictface2,adj_ang);
        agg_pts = Rotate(normalize(agg_pts),0,agg_rotate(1),agg_rotate(2)) + agg_cm;
end
        %Finding if the translated aggregate touches any other aggregate or has
        %translated inside of the center aggregate using inShape.
        clusterA = alphaShape(cluster_pts(:,1:3));
        locate = sum(inShape(clusterA,agg_pts)); 

            iter = 0;
            while locate == 0 %0 means that no points are touching (none of the points are inside)
                iter = iter+1;
                
                %The translation vector defined as the direction between 
                %center points of the fixed aggregates/cluster and translating aggregate.
                tVect = (target - agg_cm)./75; 
                agg_pts = agg_pts + tVect; %Slowly approaching towards the center
                locate = sum(inShape(clusterA,agg_pts));
                    if locate~=0
                        agg_pts = agg_pts - tVect; %If it is inside, move back a step size to remain outside.
                    end
            end
end

function result = optAngle(pts1,pts2,restrictface1,restrictface2,match_ang)
%The optimal angle to rotate the aggregates. Depending on the location of
%the aggregate/cube cell number or index. Aggregates not belonging the edge
%columns will rotate according to the adjacent aggregates while aggregates
%near the edge will rotate so that they will fit inside the two surrounding
%aggregates by matching the angles of the two sides of the aggregates to
%the sides of the two adjacent aggregates
 %Inputs:
    %pts1 - coordinates of the first aggregate
    %pts2 - coordinates of the second aggregate
    %restrictface1 - vector with the order of axis to patrition first
    %aggregate
    %strictface2 - vector with the order of axis to partrition second
    %aggregate
    %match_ang - (only for edge aggregates) rotate to match the angle of
    %the adjacent aggregates
%Outputs:
    %result - optimal angle to rotate the aggregate
a = linspace(-pi/3,pi/3,50);
optimalangle = [];
plane1 = tangentP(pts1,restrictface1,[0 0 0],false);

for ang_z = 1:length(a) %Testing each rotation to find the optimal.
    rots_z = a(ang_z);
    parfor ang_y = 1:length(a)
    rots_y = a(ang_y);
    plane2 = tangentP(pts2,restrictface2,[0, rots_y, rots_z],false);
    radt = anglebwPlanes(plane1,plane2);
    
    if ~isempty(match_ang)
        diff = abs(match_ang - radt);
        opt_ang_info = [diff, rots_y, rots_z];
    else
        opt_ang_info = [radt, rots_y, rots_z];
    end
    optimalangle = [optimalangle; opt_ang_info];
    end
end

[~,residx] = min(optimalangle(:,1));
result = optimalangle(residx,2:3);
end

function result = anglebwPlanes(normvec1,normvec2)
%Calculate the angle between two tangent plane using dot product.
%Inputs:
    %normvec1 - Normal vector of first aggregate
    %normvec2 - Normal vector of second aggregate
%Output:
    %result - the angle between the two tangent planes
magNormVec1 = norm(normvec1(1:3));
magNormVec2 = norm(normvec2(1:3));

result = acos(dot(normvec1(1:3),normvec2(1:3))./(magNormVec1.*magNormVec2));
end

function result = tangentP(pts,axis_order,angle,figurecheck)
%Calculate the normal vector of the tangent plane of the selected face of
%the aggregate. To isolate the face, the object is split by the three axes
%(x,y,z) and the face is defined by the negative or positive of the axes.
%Using those points, the extreme point is used as an "anchor" to generate
%the vectors that approximates the plane of the surface. 
%Inputs:
    %pts -  the translated aggregate coordinates with the columns
    %           representing the (x,y,z) coordinates
    %axis_order - the order of axis to partrition the aggregates
    %figurecheck - plot the aggregate with the normal plane (true - plot,
    %false - disable the figure)
%Outputs:
    %result - the coefficients of the tangent plane written as [A B C D]
    %from equation Ax+By+Cz+D=0

split1 = axis_order(1);
split2 = axis_order(2);
split3 = axis_order(3);
restrictface = axis_order(4);
pts = Rotate(normalize(pts),angle(1),angle(2),angle(3));

%Generating the tangent line of the top...
[~,idx] = sort(pts(:,split1));
sortedpts1 = pts(idx,:);
negidx = find(sortedpts1(:,split1) <= 0); 
if restrictface ==1 || restrictface == 2 
    newpts = sortedpts1(1:negidx(end,1),:); 
else
    newpts = sortedpts1(negidx(end,1)+1:end,:);
end

[~,idx] = sort(newpts(:,split2));
sortedpts2 = newpts(idx,:);
negidx = find(sortedpts2(:,split2) <= 0);
    if restrictface == 1 || restrictface == 3
    newpts = sortedpts2(1:negidx(end,1)-1,:);
    else
    newpts = sortedpts2(negidx(end,1)+1:end,:);
    end

if figurecheck == true
    plot3(newpts(:,1),newpts(:,2),newpts(:,3),'.r')
    hold on
end

%calculating the centroid
x = newpts(:,1);
y = newpts(:,2);
z = newpts(:,3);

xcm = sum(x)./length(x);
ycm = sum(y)./length(y);
zcm = sum(z)./length(z);
r0 = [xcm ycm zcm];

%Finding closest point to centroid...
order = abs(newpts(:,1:2)-r0(1:2));
[~,minidx] = min(order(:,1) + order(:,2));
r0 = newpts(minidx,:);
% plot3(r0(1),r0(2),r0(3),'.','MarkerSize',50)
% hold on

%Writing the equation for tangent plane...
%Finding the point in the same x-axis (relatively)
ptsX = [];
for i = 1:length(newpts)
    difX= newpts(i,split2) - r0(split2);
    if abs(difX) < 65
        term = newpts(i,:);
        ptsX = [ptsX; term];
    end
end
%Obtaining the vector with the most difference in Y
[~,maxY] = max(ptsX(:,split3));
[~,minY] = min(ptsX(:,split3));
WidV = ptsX(maxY,:) - ptsX(minY,:);

%Finding the point in the same y-axis (relatively)
ptsY = [];
for i = 1:length(newpts)
    difY= newpts(i,1) - r0(split3);
    if abs(difY) < 65
        term = newpts(i,:);
        ptsY = [ptsY; term];
    end
end
DepV = ptsY(end,:) - ptsY(split2,:);

Normvec = cross(DepV,WidV);
Normvec = Normvec./norm(Normvec);

%Limit x and y to only the max width/length
limxy = abs([max(newpts(:,1)), max(newpts(:,2)), min(newpts(:,1)), min(newpts(:,2))]);
limxy = max(limxy);

[X,Y] = meshgrid(-limxy:2:limxy);
%move r0 so that the z is now zmin...
if restrictface == 1 || restrictface == 2
    [~,minZ] = min(newpts(:,split1));
    r0 = newpts(minZ,:); %readjust to include the minimal point in plane
else
    [~,maxZ] = max(newpts(:,split1));
    r0 = newpts(maxZ,:);
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