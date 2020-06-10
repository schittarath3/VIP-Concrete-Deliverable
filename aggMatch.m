pts = normalize(stlread('ag1p3solid.stl').Points);
cnt = stlread('ag1p3solid.stl').ConnectivityList;
set(0,'DefaultFigureVisible','on')
TR = trimesh(cnt,pts(:,1),pts(:,2),pts(:,3));
TR.FaceAlpha = 1;
TR.EdgeColor = 'b';
axis equal

%selecting only the datapoints below the z-axis
[~,idx] = sort(pts(:,3));
sortedpts = pts(idx,:);
negidx = find(sortedpts(:,3) <= 0);
newpts = sortedpts(1:negidx(end,1),:);

%create the "mold" or allowable area for an aggregate to fit
%creating the box that fits the aggregate
x = pts(:,1);
y = pts(:,2);
z = pts(:,3);

xpt= [max(x); abs(min(x))];
ypt = [max(y); abs(min(y))];
lim = max([xpt;ypt]);
xyint = (lim+lim)./10;
[xg, yg] = meshgrid(-lim:xyint:lim); %#ok<BDSCI>
xg = xg(:);
yg = yg(:);

zg = ones(numel(xg),1);
elrep = 10;
xg = repmat(xg,elrep,1); %repeat 5 times
yg = repmat(yg,elrep,1); %repeat 5 times

%determining the intervals of the mesh...
zmax = max(newpts(:,3));
zmin = min(z);
int = (zmax - zmin) ./ (elrep -1);
zg = zg*(zmin:int:zmax); %#ok<BDSCI>
zg = zg(:);

T = delaunay(xg,yg,zg);
trimesh(T,xg,yg,zg);

figure
boxShp = alphaShape(xg,yg,zg);
in = inShape(aggShp,xg,yg,zg);
figure
xg = [xg(~in); newpts(:,1)];
yg = [yg(~in); newpts(:,2)];
zg = [zg(~in); newpts(:,3)];
alp = alphaShape(xg,yg,zg);
plot(alp)

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