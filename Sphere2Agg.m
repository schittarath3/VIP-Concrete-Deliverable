function sphAggRepo= Sphere2Agg(addAggRepo, sphereCellsrepo)
%Fill the sphere cells with the scaled aggregates from the generated
%repository containing aggregates following the GSD.
%Inputs:
%   addAggRepo - the repoistory with scaled aggregates
%   sphereCellsrepo - the repository with the sphere cells and the center
%   coordinates
%Outputs:
%   sphAggRepo - repository containing the new coordinates for the scaled
%   coordinates centered at the located sphere cells.

fields_agg = fieldnames(addAggRepo);
%Obtaining the threshold sieve size for each of the aggregates
binSz = cell2mat(sphereCellsrepo(:,2)); 

%Generating the repository of the aggregates according to their bin
%(separated by the grain/sieve size)
for aggs = 1:length(fields_agg)
    aggName = addAggRepo.(fields_agg{aggs});
    sieve_size = aggName.Diameter;
    
    fit = find(binSz>=sieve_size);
        if ~isempty(fit)
        sphereCell_fit = fit(1);
        aggRepoSort{sphereCell_fit,1}.(fields_agg{aggs}).Index = aggs;
        aggRepoSort{sphereCell_fit,1}.(fields_agg{aggs}).Diameter = aggName.Diameter;
        end
end

%Adding each of the aggregates in the list to the appropriate sphere cell
%by bins and generating the final repostitory...
for bins = 1:length(binSz)
    binSz_fields = fieldnames(aggRepoSort{bins,1});
    len_aggs = length(binSz_fields);
    
    for insertsph = 1:len_aggs
        sph_cm = cell2mat(sphereCellsrepo{bins,1}(insertsph,2));
        ang = linspace(-pi/8,pi/8,5);
        
        agg_index = aggRepoSort{bins,1}.(char(binSz_fields(aggsidxm(insertsph)))).Index;
        
        %Rewriting the coordinates for the aggregates to substitute...
        pts = addAggRepo.(fields_agg{agg_index}).OriginalPoints;
        orient = addAggRepo.(fields_agg{agg_index}).Orientation;
        pts = translate2Center(Rotate(pts,ang(orient(1)),ang(orient(2)),ang(orient(3))),sph_cm);
        
        %New repository
        newfieldname = strcat('bin',num2str(bins),'_','agg',num2str(insertsph));
        sphAggRepo.(newfieldname).Original = addAggRepo.(fields_agg{agg_index}).Original;
        sphAggRepo.(newfieldname).OriginalPoints = pts;
        sphAggRepo.(newfieldname).OriginalFaces = addAggRepo.(fields_agg{agg_index}).OriginalFaces;
        sphAggRepo.(newfieldname).Points = addAggRepo.(fields_agg{agg_index}).Points;
        sphAggRepo.(newfieldname).Faces = addAggRepo.(fields_agg{agg_index}).Faces;
        sphAggRepo.(newfieldname).Orientation = addAggRepo.(fields_agg{agg_index}).Orientation;
        sphAggRepo.(newfieldname).Diameter = addAggRepo.(fields_agg{agg_index}).Diameter;
        sphAggRepo.(newfieldname).bin = addAggRepo.(fields_agg{agg_index}).bin;
    end
end
end

function datapointsn = translate2Center(datapoints,center)
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
dcm = center - centroid; %Distance from centroid to center
for vertice = 1:length(datapoints)
    datapointsn(vertice,:) = datapoints(vertice,:) + dcm;
end
end

function nv = Rotate(pts,tx,ty,tz) 
%Rotational matrix
rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
rz = [cos(tz) -sin(tz) 0; sin(tz) cos(tz) 0; 0 0 1];
rotm = rx*ry*rz;

nv = rotm*pts';
nv = nv';
end