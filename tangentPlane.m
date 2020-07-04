%% PACKING THE AGGREGATES
clc
clear
repos = load('repos.mat');
repos = repos.repos;
fields = fieldnames(repos.myRepo);

%indexing each of the aggregates
order = ones(3,3);
iter = 0; multiplier = length(fields).^(1/3);
for j = 1:multiplier
    for i = 1:multiplier
        for k = 1:multiplier
            iter = iter + 1;
            order(k,i,j) = iter;
        end
    end
    order(:,:,j) = rot90(order(:,:,j),2);
end

figure(1) %BEFORE FIGURE
for i = 1:length(fields)
    pts = repos.myRepo.(fields{i}).Points;
    cnt = repos.myRepo.(fields{i}).Faces;
    trimesh(cnt,pts(:,1),pts(:,2),pts(:,3))
    hold on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
hold off

figure(2) %AFTER FIGURE
[pts14,cm14] = normalize(repos.myRepo.(fields{14}).Points);
pts14 = pts14 + cm14; aggpts = pts14;
cnt14 = repos.myRepo.(fields{14}).Faces;
trimesh(cnt14,pts14(:,1),pts14(:,2),pts14(:,3))
hold on

sequence = randperm(27); %auto generate a packed version of the aggregates
for idx = 1:sequence
    if idx~=14
    [pts,cm] = normalize(repos.myRepo.(fields{idx}).Points); pts = pts+cm;
    cnt = repos.myRepo.(fields{idx}).Faces;
    pts = translate2Pt(idx,pts,aggpts,cm,cm14);
    aggpts = [aggpts; pts];

    trimesh(cnt,pts(:,1),pts(:,2),pts(:,3))
    hold on
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal;
    end
end
hold off

function pts = translate2Pt(cubeidx,cubepts,aggpts,cubecm,centercm)
switch cubeidx
    case {2, 11, 20}
        restrictface1 = 'front'; rotidx = 2;
        if cubeidx == 2
            restrictface2 = 'right';
        elseif cubeidx == 11
            restrictface2 = 'back';
        elseif cubeidx == 20
            restrictface2 = 'left';
        end
        
    case {8, 17, 26}
        restrictface1 = 'back'; rotidx = 2;
        if cubeidx == 8
            restrictface2 = 'right';
        elseif cubeidx == 17
            restrictface2 = 'front';
        elseif cubeidx == 26
            restrictface2 = 'left';
        end
        
    case 13
        restrictface1 = 'bottom';
        restrictface2 = 'top';
        rotidx = 3;
        
    case 15
        restrictface1 = 'top';
        restrictface2 = 'bottom';
        rotidx = 3;
        
    case {1, 4, 7, 10, 16, 19, 22, 25}
        restrictface1 = 'bottom'; rotidx = 2;
        if cubeidx == 1 || cubeidx == 4 || cubeidx == 7
            restrictface2 = 'right';
        elseif cubeidx == 19 || cubeidx == 22 || cubeidx == 25
            restrictface2 = 'left';
        elseif cubeidx == 16
            restrictface2 = 'front';
        elseif cubeidx == 10
            restrictface2 = 'back';
        end
        
    case {3, 6, 9, 12, 18, 21, 24, 27}
        restrictface1 = 'top'; rotidx = 2;
        if cubeidx == 3 || cubeidx == 6 || cubeidx == 9
            restrictface2 = 'right';
        elseif cubeidx == 12 || cubeidx == 18 || cubeidx == 21
            restrictface2 = 'left';
        elseif cubeidx == 24
            restrictface2 = 'front';
        elseif cubeidx == 27
            restrictface2 = 'back';
        end
        
    case 23 
        restrictface1 = 'right';
        restrictface2 = 'left';
        rotidx = 2;
    case 5 
        restrictface1 = 'left';
        restrictface2 = 'right';
        rotidx = 2;
end

result = optAngle(normalize(aggpts),normalize(cubepts),restrictface1,restrictface2,rotidx);
%rotate the aggregate and return it to the orginal center
if rotidx == 2
    pts = Rotate(normalize(cubepts),0,result,0) + cubecm;
else
    pts = Rotate(normalize(cubepts),0,0,result) + cubecm;
end

%translate to the fixed point
agg1 = alphaShape(aggpts);
locate = sum(inShape(agg1,cubepts));
    while locate == 0
    [~,cubecm] = normalize(pts);
    tVect = (centercm-cubecm)./35;
    pts = pts + tVect;
    locate = sum(inShape(agg1,pts));
    
    if locate ~=0
        pts = pts - tVect; %revert back to the original location
    end
    end
end

% %volume check
% TR = triangulation(aggCubeC,aggCube);
% stlwrite(TR,'packedagg.stl') %creating an stl file
% 
% model = createpde;
% importGeometry(model,'packedagg.stl');
% mesh = generateMesh(model);
% mv = volume(mesh);
% 
% pts = stlread('packedagg.stl').Points;
% %volume of box
% Vol = (max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2)))*(max(pts(:,3)) - min(pts(:,3)));
% VolFract = mv/Vol;

%% FUNCTION CODE
function pts2 = plotaggregate(pts2,pts1,cnt2,rotangle,NormVec,option,figurecheck)
if ~isempty(rotangle)
pts2 = Rotate(pts2,angle(1),angle(2),angle(3));
end

if ~isempty(NormVec)
pts2 = Translate(pts2,pts1,option);
end

if figurecheck == true
trimesh(cnt2,pts2(:,1),pts2(:,2),pts2(:,3),'FaceAlpha',1,'EdgeColor','k','FaceColor','b');
hold on
end
end

function result = optAngle(pts1,pts2,restrictface1,restrictface2,opt)
%Inputs: agg1 - the stationary/fixed aggregate
%             agg2 - the aggregate that needs to be "attached"/shifted
a = linspace(-pi/6,pi/6);
optimalangle = [];
plane1 = tangentP(pts1,[],restrictface1,[0 0 0],false);

parfor ang = 1:length(a)
    if opt == 2 %y-axis
        rots = [0, a(ang), 0];
    else 
        rots = [0, 0, a(ang)];
    end
    plane2 = tangentP(pts2,[],restrictface2,rots,false);
    radt = anglebwPlanes(plane1,plane2);
    optimalangle = [optimalangle; radt];
end

[~,idx] = min(optimalangle);
result = a(idx);
end

function result = anglebwPlanes(normvec1,normvec2)
magNormVec1 = norm(normvec1(1:3));
magNormVec2 = norm(normvec2(1:3));

result = acos(dot(normvec1(1:3),normvec2(1:3))./(magNormVec1.*magNormVec2));
end

function result = tangentP(pts1,cnt1,option,angle,figurecheck)
%restrictface 1,2 - negative, 3,4 - positive
%restrictface 1,3 - negative to 0, 2,4 - 0 to positive

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
elseif contains(option,'front')
    split1 = 2;
    split2 = 3;
    split3 = 1;
    restrictface = 4;
elseif contains (option, 'back')
    split1 = 2;
    split2 = 3;
    split3 = 1;
    restrictface = 1;
end
pts1 = Rotate(pts1,angle(1),angle(2),angle(3));

%generating the tangent line of the top...
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

%finding closest point to centroid...
order = abs(newpts1(:,1:2)-r0(1:2));
[~,minidx] = min(order(:,1) + order(:,2));
r0 = newpts1(minidx,:);
% plot3(r0(1),r0(2),r0(3),'.','MarkerSize',50)
% hold on

%writing the equation for tangent plane...
%finding the point in the same x-axis (relatively)
ptsX = [];
for i = 1:length(newpts1)
    difX= newpts1(i,split2) - r0(split2);
    if abs(difX) < 25
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
    if abs(difY) < 25
        term = newpts1(i,:);
        ptsY = [ptsY; term];
    end
end
DepV = ptsY(end,:) - ptsY(split2,:);

Normvec = cross(DepV,WidV);
Normvec = Normvec./norm(Normvec);

%limit x and y to only the max width/length
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
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end