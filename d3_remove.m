clc
clear

im1 = imread("D:\VIP Concrete\Example\3. Binary Images\bottom_left\bl_top0076_post.png");
im2 = imread("D:\VIP Concrete\Example\3. Binary Images\bottom_left\bl_top0077_post.png");
im3 = imread("D:\VIP Concrete\Example\3. Binary Images\bottom_left\bl_top0078_post.png");
%create psuedo binary image
% mat1 = [1 1 0 0; 1 1 0 0; 0 0 0 0; 0 0 0 0];
% mat2 = [0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 1 1];
% mat3 = [0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 1 1];

%create 3D matrix
mat3D = cat(3, im1, im2, im3);

CC = bwconncomp(mat3D, 18);

pixels = CC.PixelIdxList;

a = size(CC.PixelIdxList);
num_clusters = a(2);

for cluster = 1:num_clusters
    if cluster == 36
        continue
    end
    for i = 1:length(pixels{1, cluster})
        mat3D(pixels{1, cluster}(i)) = 0;
    end
end

