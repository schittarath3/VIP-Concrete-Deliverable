clc
clear

im1 = imread("Example\3. Binary Images\bottom_left\bl_top0076_post.png");
image_size = size(im1);

mat3D = [];
base = "D:\VIP Concrete Git\Example\3. Binary Images\bottom_left\bl_top00";
for i = 54:94
    fn = strcat(base, num2str(i), '_post.png');
    a = imread(fn);
    mat3D(:,:,i-53) = a;
end

CC = bwconncomp(mat3D, 18);

pixels = CC.PixelIdxList;
a = size(CC.PixelIdxList);
num_clusters = a(2);

new_image_stack = zeros(image_size(1), image_size(2), 41);
for cluster = 1:num_clusters
    for i = 1:length(pixels{1, cluster})
        new_image_stack(pixels{1, cluster}(i)) = 1;
    end
end

image_base = "D:\VIP Concrete Analysis\Cluster Test\layer"
for i = 1:41
    img = new_image_stack(:,:,i)
    imwrite(img, strcat(image_base, num2str(i), '.png'))
end

