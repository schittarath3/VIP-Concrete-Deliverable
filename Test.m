clc
clear

ag1 = stlread('ag1.stl');
ag1_pts = ag1.Points;
ag1_cl = ag1.ConnectivityList;

%% Generating Repository of Aggregates
%Load in all of the aggregates and shrink them by percentage (set this as a
%variable so we can experiment later) - to shrink, use the scale matrix..

%Obtain the volume of the aggregates mesh (volume that has been shrunk) and
%the volume of the box by finding the min,max values, store the width,
%length, and height values...

%Apply a set of transformations to all of the aggregates and do the same
%thing (store the values of the volume mesh and the width, length, and
%height)
    %Transformations: Limit the angles? [-pi/12,pi/12], set the rotation
    %matrices and apply a for loop to go through all of them (rip)...
        %for theta_x = -pi/12: pi/12
            %for theta_y = -pi/12: pi/12
                %for theta_z = -pi/12 : pi/12...

%OUTPUT: A matrix with all of the aggregates ID, include properties: Mesh
%Volume, Box Volume (for packing), Width, Length, and Height (dim of box
%volume)

%% FUNCTION 1: Normalize Each Aggregate
%Obtain the center of each aggregate
x = ag1_pts(:,1);
y = ag1_pts(:,2);
z = ag1_pts(:,3);

xcm = sum(x)./length(x);
ycm = sum(y)./length(y);
zcm = sum(z)./length(z);
centroid = [xcm ycm zcm];

%Obtain the matrix with the distance of each vertices to the center
ag1_dpts = ag1_pts;
dcm = centroid - [0 0 0]; %Distance from centroid to (0,0,0)
for vertice = 1:length(ag1_pts)
    ag1_dpts(vertice,:) = [ag1_pts(vertice,:) - dcm];
end

%Obtain the matrix with the new coordinates now center at (0,0,0)
figure
trimesh(ag1_cl, ag1_dpts(:,1), ag1_dpts(:,2), ag1_dpts(:,3));

%% FUNCTION 2: Generating the Volume of Each Aggregate
width = abs(max(ag1_dpts(:,1)) - min(ag1_dpts(:,1)));
length = abs(max(ag1_dpts(:,2)) - min(ag1_dpts(:,2)));
height = abs(max(ag1_dpts(:,3)) - min(ag1_dpts(:,3));
V = width.*length.*height;

%% PACKING THE AGGREGATES
%Feed the algorithm with the first 10-20 aggregates that MUST be contained
%inside. Start with one aggregate with center at (0,0,0), then keep adding
%all of the other aggregates one-by-one. First define the min and max
%dimensions of the box then go into the repository and find the ones that
%best fit. Switch the starting aggregate (for loop here). 
    %for i = 1: # of aggregates MUST be in the box...
        %Place first aggregate inside the box (already at 0,0,0)
        
        %Find distance from the most left of aggregate to most left of the
        %container then go to repository to find the aggregate with closest
        %dimension that fits...
        
        %Place second aggregate inside the box and store the new
        %coordinates (apply the shift), blah blah...
        
        %disp('Success') or ('Fail') if the aggregates do not all fit in?