function VolFract = tangentPlane(repository)
repos = load(repository);
fields = fieldnames(repos.myRepo);
angles = linspace(-pi/8,pi/8,5);

%create a matrix to index each of the aggregate...
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

%generating the columns of each aggregates
for col = 1:3
    for agg = 1:3
        ao = order(1:3,agg,col);
        ang = linspace(-pi/8,pi/8,5);
        rot1 = repos.myRepo.(fields{ao(3)}).Orientation; 
        rot2 = repos.myRepo.(fields{ao(2)}).Orientation;
        rot3 = repos.myRepo.(fields{ao(1)}).Orientation;
        ptsn1 = Rotate(normalize(stlread(strcat(repos.myRepo.(fields{ao(3)}).Original,'.stl')).Points),ang(rot1(1)),ang(rot1(2)),ang(rot1(3)));
        ptsn2 = Rotate(normalize(stlread(strcat(repos.myRepo.(fields{ao(2)}).Original,'.stl')).Points),ang(rot2(1)),ang(rot2(2)),ang(rot2(3)));
        ptsn3 = Rotate(normalize(stlread(strcat(repos.myRepo.(fields{ao(1)}).Original,'.stl')).Points),ang(rot3(1)),ang(rot3(2)),ang(rot3(3)));
        
        cntn1 = stlread(strcat(repos.myRepo.(fields{ao(3)}).Original,'.stl')).ConnectivityList;
        cntn2 = stlread(strcat(repos.myRepo.(fields{ao(2)}).Original,'.stl')).ConnectivityList;
        cntn3 = stlread(strcat(repos.myRepo.(fields{ao(1)}).Original,'.stl')).ConnectivityList;
        
        rottop = optAngle(ptsn2,ptsn1,'bottom','top',2);
            Normvec = tangentPlane(ptsn1,cntn1,'top',[0 rottop 0],false);
            newpts = [plotaggregate(ptsn1,ptsn2,cntn1,[0 rottop 0],Normvec,'top',false); ptsn2];
            
        rotbot = optAngle(ptsn2,ptsn3,'top','bottom',2);
            Normvec = tangentPlane(ptsn3,cntn3,'bottom',[0 rotbot 0],false);
            newpts = [plotaggregate(ptsn3,newpts,cntn3,[0 rotbot 0],Normvec,'bottom',false); newpts];
            
        plotaggregate(ptsn2,[],cntn2,[],[],[],true);
        
        cntpts = [cntn3; cntn1+length(ptsn3); cntn2+length([ptsn1; ptsn3])];
        aggStacks.(strcat('stack',num2str(ao(1)./3))) = newpts;
        aggStacksC.(strcat('stack',num2str(ao(1)./3))) = cntpts;
    end
end

%generating each of the layers...
for layer = 1:3
    stackidx = order(:,layer,1);
    stack1 = aggStacks.(strcat('stack',num2str(stackidx(3))));
    stackcn1 = aggStacksC.(strcat('stack',num2str(stackidx(3))));
    stack2 = aggStacks.(strcat('stack',num2str(stackidx(2))));
    stackcn2 = aggStacksC.(strcat('stack',num2str(stackidx(2))));
    stack3 = aggStacks.(strcat('stack',num2str(stackidx(1))));
    stackcn3 = aggStacksC.(strcat('stack',num2str(stackidx(1))));
    
    rotstack = optAngle(stack2,stack1,'right','left',2);
    Normvec = tangentPlane(stack1,stack1,'left',[0 0 rotstack],false);
    newstack = [plotaggregate(stack1,stack2,stackcn1,[0 0 rotstack],Normvec,'left',false); stack2];
    
    rotstack = optAngle(stack2,stack3,'left','right',2);
    Normvec = tangentPlane(stack1,stack1,'right',[0 0 rotstack],false);
    newstack = [plotaggregate(stack3,newstack,stackcn3,[0 0 rotstack],Normvec,'right',false); newstack];
    
    plotaggregate(stack2,[],stackcn2,[],[],[],true);
    
    cntpts = [stackcn3; stackcn1+length(stack3); stackcn2+length([stack1; stack3])];
    aggLayer.(strcat('layer',num2str(layer))) = newstack;
    aggLayerC.(strcat('layer',num2str(layer))) = cntpts;
end

for cube = 1:3
    layer1 = aggLayer.(strcat('layer','3')); layercn1 = aggLayerC.(strcat('layer','3'));
    layer2 = aggLayer.(strcat('layer','2')); layercn2 = aggLayerC.(strcat('layer','2'));
    layer3 = aggLayer.(strcat('layer','1')); layercn3 = aggLayerC.(strcat('layer','1'));
    figure(1)
    aggCube = [plotaggregate(layer1,layer2,layercn1,[0 0 0],Normvec,'front',true); layer2];
    hold on
    aggCube = [plotaggregate(layer3,aggCube,layercn3,[0 0 0],Normvec,'back',true);aggCube];
    hold on
    plotaggregate(layer2,[],layercn2,[],[],[],true);
    axis equal
    hold off
    
    aggCubeC = [layercn3; layercn1+length(layer3); layercn2+length([layer1; layer3])];
end

%VolumeCheck
TR = triangulation(aggCubeC,aggCube);
stlwrite(TR,'packedagg.stl')

model = createpde;
importGeometry(model,'packedagg.stl');
mesh = generateMesh(model);
mv = volume(mesh);

pts = stlread('packedagg.stl').Points;
%volume of box
Vol = (max(pts(:,1)) - min(pts(:,1)))*(max(pts(:,2)) - min(pts(:,2)))*(max(pts(:,3)) - min(pts(:,3)));
VolFract = mv/Vol;
end

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
a = linspace(-pi/3,pi/3);
optimalangle = [];
plane1 = tangentPlane(pts1,[],restrictface1,[0 0 0],false);

parfor ang = 1:length(a)
    if opt == 2 %y-axis
        rots = [0, a(ang), 0];
    else 
        rots = [0, 0, a(ang)];
    end
    plane2 = tangentPlane(pts2,[],restrictface2,rots,false);
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

function result = tangentPlane(pts1,cnt1,option,angle,figurecheck)
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

function datapointsn = Translate(datapoints1,datapoints2,option)
agg2 = alphaShape(datapoints2); %stationary
alpha = 0;
StepSz = 1;
datapointsn = datapoints1;
locate = sum(inShape(agg2,datapoints1));

while locate~=0
if contains(option,'left')
    datapointsn(:,1) = datapointsn(:,1) + -1*(alpha+StepSz);
elseif contains(option,'right')
    datapointsn(:,1) = datapointsn(:,1) + 1*(alpha+StepSz);
elseif contains(option,'bottom')
    datapointsn(:,3) = datapointsn(:,3) + 1*(alpha+StepSz);
elseif contains(option,'top')
    datapointsn(:,3) = datapointsn(:,3) + -1*(alpha+StepSz);
elseif contains(option,'front')
    datapointsn(:,2) = datapointsn(:,2) + 1*(alpha+StepSz);
elseif contains(option,'back')
    datapointsn(:,2) = datapointsn(:,2) + -1*(alpha+StepSz);
end
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
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end