clc 
clear

filename = 'ag1p.stl';
agg = stlread(filename);

%generating the points and connectivity
points = normalize(agg.Points);

%finding volume of mesh
model = createpde;
importGeometry(model,filename);
mesh = generateMesh(model);
Vmesh = volume(mesh);

%finding original fraction
Vbox = getVolBox(points,0);
oldVfract = Vmesh/Vbox;

%simulated annealing optimization
t = sym('t', [3 1]);
Rx = symfun([1 0 0; 0 cos(t(1)) -sin(t(1)); 0 sin(t(1)) cos(t(1))],t);
Ry = symfun([cos(t(2)) 0 sin(t(2)); 0 1 0; -sin(t(2)) 0 cos(t(2))],t);
Rz = symfun([cos(t(3)) -sin(t(3)) 0; sin(t(3)) cos(t(3)) 0; 0 0 1],t);
R = Rx*Ry*Rz;
pnew = R*points';

theta = [-pi/12; -pi/12; 0];
pnew2 = subs(pnew,{t(1),t(2),t(3)},{theta(1),theta(2),theta(3)});
pnew2 = formula(pnew2);

%finding the index of the max and min (if plugged in said points)...
[mxw, idx_mxw] = max(pnew2(1,:));
[mnw, idx_mnw] = min(pnew2(1,:));
[mxl, idx_mxl] = max(pnew2(2,:));
[mnl, idx_mnl] = min(pnew2(2,:));
[mxh, idx_mxh] = max(pnew2(3,:));
[mnh, idx_mnh] = min(pnew2(3,:));

%reverse to obtain new measurement...
pnew = formula(pnew);
w = pnew(1,idx_mxw) - pnew(1,idx_mnw);
l = pnew(2,idx_mxl) - pnew(2,idx_mnl);
h = pnew(3,idx_mxh) - pnew(3,idx_mnh);
V = matlabFunction(formula(w*l*h));
Volume = @(t) V(t(1),t(2),t(3));

%using simulated annealing to find optimal angle...
x0 = [0, 0, 0];
lb = [-pi/12, -pi/12, -pi/12];
ub = [pi/12, pi/12, pi/12];
result = simulannealbnd(Volume,x0,lb,ub);

tx = result(1);
ty = result(2);
tz = result(3);
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
Ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
R = Rx*Ry*Rz;
dtpoints = R*points';

box = getVolBox(dtpoints,0);
newVfract = Vmesh/box;

%% Functions
function [Vbox, l, w, h] = getVolBox(datapoints,opt)
%Input: datapoints - written as column vectors [x y z];
%Output: V - the volume of the container/box
%            l, w, h - dimensions of the box (length, width, height)
w = max(datapoints(:,1)) - min(datapoints(:,1));
l = max(datapoints(:,2)) - min(datapoints(:,2));
h = max(datapoints(:,3)) - min(datapoints(:,3));
Vbox = w.*l.*h;

xpt = [min(datapoints(:,1)); max(datapoints(:,1))];
ypt = [min(datapoints(:,2)); max(datapoints(:,2))];
zpt = [min(datapoints(:,3)); max(datapoints(:,3))];
    
if opt == 1
    box_ptstop = [xpt(1) ypt(2) zpt(2); xpt(1) ypt(1) zpt(2); xpt(2) ypt(2) zpt(2); xpt(2) ypt(1) zpt(2)];
    box_ptsbot = [xpt(1) ypt(2) zpt(1); xpt(1) ypt(1) zpt(1); xpt(2) ypt(2) zpt(1); xpt(2) ypt(1) zpt(1)];
    box_pts = [box_ptstop; box_ptsbot];
    box_connect = [1 2 3; 2 3 4; 5 6 7; 6 7 8; 2 4 6; 6 4 8; 4 3 8; 8 3 7; 1 3 5; 3 5 7; 1 2 5; 2 5 6];
    trimesh(box_connect,box_pts(:,1),box_pts(:,2),box_pts(:,3),'FaceAlpha',0.5)
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