%find min and maxes and make a new bounding box
function volumeFrac = volFracCheck(aggRepo)
%Gets volume fraction of aggregate within the cube
%Inputs: 
%   aggRepo: Struct of aggregates of form aggName ->cubeNum,Faces,Points
%   cubeLength: Length of cube
%Outputs:
%   Volume fraction of aggregates within cube as float value
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    totalAggVol = 0;
    minX = 100000;
    maxX = -100000;
    minY = 100000;
    maxY = -100000;
    minZ = 100000;
    maxZ = -100000;
    for i = 1:numAgg %calculate volume for each aggregate and totals
        p = aggRepo.(aggNames{i}).Points';
       
        curMinX = min([p(:,1)]);
        if curMinX < minX
            minX = curMinX;
        end
        
        curMaxX = max([p(:,1)]);
        if curMaxX > maxX
            maxX = curMaxX;
        end
        
        curMinY = min([p(:,2)]);
        if curMinY < minY
            minY = curMinY;
        end
        
        curMaxY = max([p(:,2)]);
        if curMaxY > maxY
            maxY = curMaxY;
        end
        
        curMinZ = min([p(:,1)]);
        if curMinZ < minZ
            minZ = curMinZ;
        end
        
        curMaxZ = max([p(:,1)]);
        if curMaxZ > maxZ
            maxZ = curMaxZ;
        end
        
        f = aggRepo.(aggNames{i}).Faces';
        [volume, area] = stlVolume(p,f);
        totalAggVol = totalAggVol + volume
    end
    
    xLen = maxX - minX;
    yLen = maxY - minY;
    zLen = maxZ - minZ;
    
    boundingVol = xLen * yLen * zLen;
    
    volumeFrac = -totalAggVol/boundingVol;
end

function [totalVolume,totalArea] = stlVolume(p,t)
% Given a surface triangulation, compute the volume enclosed using
% divergence theorem.
% Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
% Input: p: (3xnPoints), t: (3xnTriangles)
% Output: total volume enclosed, and total area of surface  
% Author: K. Suresh; suresh@engr.wisc.edu

% Compute the vectors d13 and d12
d13= [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
cr = cross(d13,d12,1);%cross-product (vectorized)
area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);% Area of each triangle
totalArea = sum(area);
crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
nz = -cr(3,:)./crNorm;% z component of normal for each triangle
volume = area.*zMean.*nz; % contribution of each triangle
totalVolume = sum(volume);%divergence theorem
end