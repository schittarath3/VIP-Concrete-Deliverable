cnt = repos.ag11p3.Faces;
a = linspace(-pi/8,pi/8,5);
ang = [1,1,1];
% opts = aggRepo.ag31p3_525.OriginalPoints;
pts = Rotate(repos.ag11p3.Vertices,a(ang(1)),a(ang(2)),a(ang(3)));
opts = Rotate(repos.ag11p3.OriginalPoints,a(ang(1)),a(ang(2)),a(ang(3)));
ocnt = repos.ag11p3.OriginalFaces;
pts1 = repos.ag11p3.Orientation.tx_indx1ty_indx1tz_indx1;

% cnt = aggRepo.ag16Bp3_121.OriginalFaces;
% pts = aggRepo.ag16Bp3_121.OriginalPoints;
% ocnt = aggRepo.ag16Bp3_121.Faces;
% opts = aggRepo.ag16Bp3_121.Points;

figure(1)
% trimesh(ocnt,opts(:,1),opts(:,2),opts(:,3),'EdgeColor','blue');
% hold on
trimesh(cnt,pts(:,1),pts(:,2),pts(:,3),'FaceAlpha',0,'EdgeColor','red')
hold on
trimesh(cnt,pts1(:,1),pts1(:,2),pts1(:,3))
axis equal
hold off

function nv = Rotate(pts,tx,ty,tz) 
%Rotational matrix
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
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