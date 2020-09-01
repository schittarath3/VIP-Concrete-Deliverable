tic
[x1,y1,z1] = sphere(75);
x1 = x1(:);
y1 = y1(:);
z1 = z1(:);
P = [x1 y1 z1];
P = unique(P,'rows');

[x2,y2,z2] = sphere(125);
x2 = x2(:)+2;
y2 = y2(:)+1;
z2 = z2(:)+1;
Q = [x2 y2 z2];
Q = unique(Q,'rows');

% shp1 = alphaShape(P(:,1),P(:,2),P(:,3),1);
% shp2 = alphaShape(Q(:,1),Q(:,2),Q(:,3),1);
% figure(1)
% plot(shp1)
% hold on
% plot(shp2)
% axis equal

[~,cm] = normalize(Q,[0,0,0]);
tf = sum(inShape(shp1,Q(:,1),Q(:,2),Q(:,3)));
iter = 0;
stepSz = 1;
while tf == 0
    iter = iter+1;
    newQ = Scale(Q,stepSz,stepSz,stepSz);
    newQ = normalize(newQ,cm);
     tf = sum(inShape(shp1,newQ(:,1),newQ(:,2),newQ(:,3)));
    
    if tf ~=0
        break
    elseif iter == 25
        break
    elseif tf == 0
        Q = newQ;
        stepSz = stepSz + .0025;
        
%         figure(1)
%         shp1 = alphaShape(P(:,1),P(:,2),P(:,3),1);
%         shp2 = alphaShape(Q(:,1),Q(:,2),Q(:,3),1);
%         plot(shp1,'FaceAlpha',0)
%         hold on
%         plot(shp2,'FaceAlpha',0,'EdgeColor','red')
%         axis equal
%         hold off
    end
end
toc

function nv = Scale(pts,sx,sy,sz) 
%Rotational matrix
scalem = diag([sx sy sz]);

nv = scalem*pts';
nv = nv';
end

function [datapointsn, centroid] = normalize(datapoints,target)
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
dcm = target  - centroid; %Distance from centroid to (0,0,0)
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = datapoints(vertice,:) + dcm;
end
end