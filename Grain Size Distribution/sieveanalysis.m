%Sieve analysis/particle distribution analysis of the 2D aggregates using
%the generated binary images. The calculations can be referenced to the
%research paper "Image Analysis Techniques on Evaluation of Particle Size 
%Distribution of Gravel by G.H.A. Kumara, et al.

%sieveanalysis.m must be placed directly inside the same folder with the
%concrete scan images to read correctly. 
clc; 
clear;
files = dir();
dirFlags = [files.isdir];
subFolders = files(dirFlags);

%% 
%Empty matrices to store information...
results = []; 
totpassInfo = []; 
totVolume = []; 
grainSzm = [];
    
%Workflow - read each of the images, compile into one data and plot the
%grain size distribution in figure.
for folder = 1:length(subFolders)
    try
    imageFold =  subFolders(folder).name;
    fileList = dir(strcat(imageFold,'/*.png'));
    
    scale = .10769;
    numSieve = 40;
    sieveSz = linspace(0,45/scale,numSieve); %Size of sieve in equal intervals
    
    for files = 1:length(fileList)
        filename = strcat(fileList(files).folder,'\',fileList(files).name);
        [passInfo, totVol, grsize] = processImg(filename,sieveSz);
        grainSzm = [grainSzm, grsize'];
        totpassInfo = [totpassInfo, passInfo];
        totVolume = [totVolume, totVol];
    end

    catch
        disp('not a valid folder')
    end
end

%% 
%Figures
%GSD
figure(1)
results = (sum(totpassInfo,2)./sum(totVolume)) * 100;
semilogx(sieveSz,results,'-*')
ylim([0 100])
title('Grain Size Distribution Sieve Analysis')
xlabel('Grain Size (pixels)')
ylabel('Cummulative Passing (%), Volume')
grid on

%Histogram of the Particle Sizes in Bins
figure(2)
histogram(grainSzm)
xlabel('Grain Size (pixels)')
ylabel('Grain Count')
title('Particle Size Distribution')

grainsizeresults = [results sieveSz'];

%% 
%Functions
function [pass, totVol, grainsize] = processImg(filename,sieveSz)
BW=imread(filename); 
BW=imfill(BW,'holes');
BW=wiener2(BW,[5 5]);

%Get properties of individual particles
props = regionprops('table',BW,'Centroid','MajorAxisLength',...
'MinorAxisLength','ConvexArea','Area','Perimeter');

%scale = .10769;
scale = 1;
major=table2array(props(:,{'MajorAxisLength'}))*scale;
minor=table2array(props(:,{'MinorAxisLength'}))*scale;
area=table2array(props(:,{'Area'}))*scale^2;

grainsize = sqrt(0.5*(minor.^2 + major.^2));
volumeest = (4/3).*area.*(1/2).*(0.6837*minor);

%Grain size vs. percent passing
pass = zeros(length(sieveSz),length(grainsize)); %matrix storing success in passing
for Sieve = 1:length(sieveSz) %starting with first sieve size
    for grain = 1:length(grainsize) 
        if grainsize(grain) < sieveSz(Sieve)
            pass(Sieve,grain) = volumeest(grain);
        end
    end
end
totVol = sum(volumeest);
end