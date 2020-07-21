clc;clear;

files = dir('*.stl');

angles = linspace(0,2*pi,20);
grainSz = zeros(length(angles),length(angles));


for num = 1:length(files)
for angY = 1:length(angles)
for angX = 1:length(angles)
    filename = files(num).name;
    
    if angY == 1 && angX == 1
    %volume
    model = createpde;
    importGeometry(model,filename);
    mesh = generateMesh(model);
    repos.(filename(1:end-4)).Volume  = volume(mesh);
    end
    
    pts = Rotate(normalize(stlread(filename).Points),angles(angX),angles(angY),0);
%     cnt = stlread(files(1).name).ConnectivityList;
%     
%     figure(1)
%     trimesh(cnt,pts(:,1),pts(:,2),pts(:,3))
%     axis equal
%     hold on
    [Xlength1,maxX] = max(pts(:,1));
    [Xlength2,minX] = min(pts(:,1));
    [Ylength1,maxY] = max(pts(:,2));
    [Ylength2,minY] = min(pts(:,2));

    grainSz(angX,angY) = sqrt(0.5*((Xlength1-Xlength2).^2 + (Ylength1-Ylength2).^2));
%     view(0,90)
%     plot3(pts(maxX,1),pts(maxX,2),pts(maxX,3),'r*');
%     plot3(pts(maxY,1),pts(maxY,2),pts(maxY,3),'b*');
%     plot3(pts(minY,1),pts(minY,2),pts(minY,3),'g*');
%     plot3(pts(minX,1),pts(minX,2),pts(minX,3),'m*');
%     xlabel('x')
%     ylabel('y')
%     hold off
end
end
repos.(filename(1:end-4)).GrainSize = min(grainSz,[],'all');
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