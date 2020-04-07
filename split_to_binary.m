function results = myimfcn(im)
%Image Processing Function
%
% IM      - Input image.
% RESULTS - A scalar structure with the processing results.
%

%--------------------------------------------------------------------------
% Auto-generated by imageBatchProcessor App. 
%
% When used by the App, this function will be called for every input image
% file automatically. IM contains the input image as a matrix. RESULTS is a
% scalar structure containing the results of this processing function.
%
%--------------------------------------------------------------------------



% Replace the sample below with your code----------------------------------

im = im(:,:,1); 

%wiener denoise
img2 = wiener2(im, [4 4]);

%apply CLAHE and gaussian blur so that binary image will be sharper
img3 = adapthisteq(img2);
img4 = imgaussfilt(img3, 2);

%create disk-shaped structuring element, then subtract from filtered image
%to remove illumination effects
diskSize = 75; 
SE = strel('disk',diskSize); 
nhood = ~(SE.Neighborhood); 
se = strel(nhood); 
bg_img = imopen(img4, se); 
background = imgaussfilt(bg_img, 3);
co_img = img4 - bg_img; 

%create binary image
bw = imbinarize(co_img); 

%remove artifact objects and fill in holes in binary image
pixelAreaToRemove = 600; 
bw = bwareaopen(bw, pixelAreaToRemove); 
bw = ~bwareaopen(~bw, 500);

%find connected components for separation
D = bwdist(~bw);
D = imcomplement(D);
l = watershed(D);
bw2 = bw;
bw2(l == 0) = 0;
mask = imextendedmin(D, 3);
D2 = imimposemin(D, mask);
L2 = watershed(D2);
bw3 = bw;
bw3(L2 == 0) = 0;
bw4 = imerode(bw3, strel('disk', 1))


results.post     = bw4;

%--------------------------------------------------------------------------
