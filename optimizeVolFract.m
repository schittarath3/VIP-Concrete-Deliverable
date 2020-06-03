function [newpoints, perchange] = optimizeVolFract(filename)
agg = stlread(filename);

%generating the points and connectivity
pts = normalize(agg.Points);

%finding volume of mesh
model = createpde;
importGeometry(model,filename);
mesh = generateMesh(model);
Vmesh = volume(mesh);

%finding original fraction
Vbox = getVolBox(points,0);
oldVfract = Vmesh/Vbox;

%defining the initial guesses...
tol = 10^-6;
t0 = .01;
delta = 50;
alp = 1/8;
t = [-pi/6; -pi/12; 0];

%defining the volume...
[V, ptmax, ptmin] = Volume(pts,t);
dV = dVolume(ptmax,ptmin,t);

for iter = 1:1e4
    tnew = ones(3,1);
    for i = 1:3 %applying random transformation
        R = randn;
        tnew(i,1) = t(i,1) - alp.*dV + t0.*R.*delta;
    end
    [Vnew, ptmaxnew, ptminnew] = Volume(pts,tnew);
    rho = min(exp((V-Vnew)./t0),1); %acceptance rate
    rn = rand; 

    if Vnew < V %automatically accept best solution
        t = tnew;
        V = Vnew;
        dV = dVolume(ptmaxnew,ptminnew,t);
    else %deciding whether or not accept worse solution
        if rn < rho
            t = tnew;
            [V, ptmaxnew, ptminnew,newpoints] = Volume(pts,t);
            dV = dVolume(ptmaxnew,ptminnew,t);
        end
    end
end
newVFract = Vmesh/Volume(pts,t);
perchange = newVFract/oldVFract * 100;
end

%% Functions

%derivative of volume...
function dVol = dVolume(ptmax,ptmin,theta)
t = sym('t', [3 1]);
Rx = symfun([1 0 0; 0 cos(t(1)) -sin(t(1)); 0 sin(t(1)) cos(t(1))],t);
Ry = symfun([cos(t(2)) 0 sin(t(2)); 0 1 0; -sin(t(2)) 0 cos(t(2))],t);
Rz = symfun([cos(t(3)) -sin(t(3)) 0; sin(t(3)) cos(t(3)) 0; 0 0 1],t);
R = Rx*Ry*Rz;

pmin = sym('pmin', [3 1]);
pnmin = formula(symfun(R*[pmin(1); pmin(2); pmin(3)],[pmin; t]));
xmin = sum(pnmin(1,:));
ymin = sum(pnmin(2,:));
zmin = sum(pnmin(3,:));

pmax = sym('pmax', [3 1]);
pnmax = formula(symfun(R*[pmax(1); pmax(2); pmax(3)],[pmax; t]));
xmax = sum(pnmax(1,:));
ymax = sum(pnmax(2,:));
zmax = sum(pnmax(3,:));

V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin);

dVx = diff(V,t(1));
dVy = diff(V,t(2));
dVz = diff(V,t(3));
dV = norm([dVx; dVy; dVz]);

dVCal = matlabFunction(dV);
dVCal = @(pmax,pmin,t) dVCal(pmax(1),pmax(2),pmax(3),pmin(1),pmin(2),pmin(3),t(1),t(2),t(3));

dVol = dVCal(ptmax,ptmin,theta);
end 

%volume calculation
function [V, ptmax, ptmin, newpts] = Volume(pts,t)
Rx = [1 0 0; 0 cos(t(1)) -sin(t(1)); 0 sin(t(1)) cos(t(1))];
Ry = [cos(t(2)) 0 sin(t(2)); 0 1 0; -sin(t(2)) 0 cos(t(2))];
Rz = [cos(t(3)) -sin(t(3)) 0; sin(t(3)) cos(t(3)) 0; 0 0 1];
R = Rx*Ry*Rz;

newpts = R*pts';

ptmax = [max(newpts(1,:)); max(newpts(2,:)); max(newpts(3,:))];
ptmin = [min(newpts(1,:)); min(newpts(2,:)); min(newpts(3,:))];

w = ptmax(1) - ptmin(1);
l = ptmax(2) - ptmin(2);
h = ptmax(3) - ptmin(3);

V = w*l*h; %function
end

%normalize points
function datapointsn = normalize(datapoints)
%obtain the center of each aggregate
x = datapoints(:,1);
y = datapoints(:,2);
z = datapoints(:,3);

xcm = sum(x)./length(x);
ycm = sum(y)./length(y);
zcm = sum(z)./length(z);
centroid = [xcm ycm zcm];

%obtain the matrix with the distance of each vertices to the center
datapointsn = datapoints;
dcm = [0, 0, 0] - centroid; %Distance from centroid to (0,0,0)
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = datapoints(vertice,:) + dcm;
end
end