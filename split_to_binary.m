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
D = bwdist(~bw); %compute distance transformation
D = imcomplement(D);
l = watershed(D); %watershed transformation
bw2 = bw;
bw2(l == 0) = 0;
%binary image is segmented at every minimum (oversegmented)
%so we use minima imposition
mask = imextendedmin(D, 3);
D2 = imimposemin(D, mask);
L2 = watershed(D2);
bw3 = bw;
bw3(L2 == 0) = 0;
%increase watershed line width from 1px to 3px
bw4 = imerode(bw3, strel('disk', 1));

%figure(1)
%title('before refinement');
%imshow(bw4);

%refining the images to remove aggregates near edges
for rot = 1:4 %rotating image for each side 
    bw4 = imrotate(bw4,90);
    aggregates = bwconncomp(bw4).PixelIdxList; %group aggregates by PIXELS
    for aggregate = 1:10 %eliminate first ten aggregates
        bw4(aggregates{aggregate}) = 0; 
    end %stop removing
end %stop rotating

%removing small aggregates/clean image...
bw5 = bw4; 
numPixels = cellfun(@numel,bwconncomp(bw5).PixelIdxList);
for aggregate = 1:length(numPixels)
    pixel_min = 150; %pixels too small
    if numPixels(aggregate) < pixel_min
        bw5(bwconncomp(bw4).PixelIdxList{aggregate}) = 0;
    end %stop eliminating
end %stop evaluation

%figure(2)
%title('after refininement')
%imshow(bw5)

results.post = bw5;
end

%--------------------------------------------------------------------------
