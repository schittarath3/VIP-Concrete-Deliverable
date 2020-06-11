function tangentPlane(filename,restrictface)
pts1 = normalize(stlread(filename).Points);

%generating the tangent line of the top...
[~,idx] = sort(pts1(:,3));
sortedpts1 = pts1(idx,:);
negidx = find(sortedpts1(:,3) <= 0); %gathering all the points below z=0
if restrictface ==1 || restrictface == 2 %below Z=0
    newpts1 = sortedpts1(1:negidx(end,1),:); %forming new pts below z=0
else
    newpts1 = sortedpts1(negidx(end,1)+1:end,:); %forming new pts above z=0
end

%finding if there are any flat faces
lowz = unique(round(newpts1(:,3),3));
lowzidx = find(newpts1(:,3) > lowz(1));

%removing flat face to generate representative plane...
%rotating the aggregate until the min x is more than 0...
flatface = newpts1(1:lowzidx-1,:);

if length(flatface) > 10 %needs at least 10 points for flat face
dTheta = 1; 
stepSz = .25;
if ~isempty(flatface) %there is flat face
    while min(flatface(:,1)) < 0
        dTheta = dTheta + stepSz;
        flatface = Rotate(flatface,0,0,-pi./dTheta);
    end
    pts1 = Rotate(pts1,0,0,-pi/dTheta);
    newpts1 = Rotate(newpts1,0,0,-pi/dTheta);
end
end
cnt1 = stlread(filename).ConnectivityList;
figure
mesh1 = trimesh(cnt1,pts1(:,1),pts1(:,2),pts1(:,3));
mesh1.FaceAlpha = 1;
mesh1.EdgeColor = 'b';
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hold on

%which tangent plane needs to be find depends on the region...
[~,idx] = sort(newpts1(:,1));
sortedpts2 = newpts1(idx,:);
negidx = find(sortedpts2(:,1) <= 0);
    if restrictface == 1 || restrictface == 3
    newpts1 = sortedpts2(1:negidx(end,1)-1,:);
    else
    newpts1 = sortedpts2(negidx(end,1)+1:end,:);
    end
plot3(newpts1(:,1),newpts1(:,2),newpts1(:,3),'.r')
hold on

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
plot3(r0(1),r0(2),r0(3),'.','MarkerSize',50)
hold on

%writing the equation for tangent plane...
%finding the point in the same x-axis (relatively)
ptsX = [];
for i = 1:length(newpts1)
    difX= newpts1(i,1) - r0(1);
    if abs(difX) < 8
        term = newpts1(i,:);
        ptsX = [ptsX; term];
    end
end
%obtaining the vector with the most difference in Y
[~,maxY] = max(ptsX(:,2));
[~,minY] = min(ptsX(:,2));
WidV = ptsX(maxY,:) - ptsX(minY,:);

%finding the point in the same y-axis (relatively)
ptsY = [];
for i = 1:length(newpts1)
    difY= newpts1(i,2) - r0(2);
    if abs(difY) < 8
        term = newpts1(i,:);
        ptsY = [ptsY; term];
    end
end
DepV = ptsY(end,:) - ptsY(1,:);

Normvec = cross(DepV,WidV);
Normvec = Normvec./norm(Normvec);
[X,Y] = meshgrid(-45:2:45);
%move r0 so that the z is now zmin...
if restrictface == 1 || restrictface == 2
    [~,minZ] = min(newpts1(:,3));
    r0 = newpts1(minZ,:); %readjust to include the minimal point in plane
else
    [~,maxZ] = max(newpts1(:,3));
    r0 = newpts1(maxZ,:);
end
D = -1*(Normvec(1)*r0(1) + Normvec(2)*r0(2) + Normvec(3)*r0(3));
Z = -1/Normvec(3)*(Normvec(1)*X + Normvec(2)*Y + D);
surf(X,Y,Z,'FaceColor','w');
hold off
axis equal
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

function nv = Rotate(pts,tx,ty,tz)
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end