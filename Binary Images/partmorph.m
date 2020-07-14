clc; clear;
files = dir();
dirFlags = [files.isdir];
subFolders = files(dirFlags);

results = []; totpassInfo = []; totVolume = [];
    
for folder = 1:length(subFolders)
    try
    imageFold =  subFolders(folder).name;
    fileList = dir(strcat(imageFold,'/*.png'));
    
    scale = .10769;
    numSieve = 10;
    sieveSz = linspace(0,45/scale,numSieve); %size of sieve in equal intervals
    
    for files = 1:length(fileList)
        filename = strcat(fileList(files).folder,'\',fileList(files).name);
        [passInfo, totVol] = processImg(filename,sieveSz);
        totpassInfo = [totpassInfo, passInfo];
        totVolume = [totVolume, totVol];
    end

    %x-intervals are the sieveSz
    %y-intervals are the percent passing: the number of aggregates that pass
    %through divided by total number of aggregates (size of grainsize)
%     results = zeros(numSieve,1);
%     for k = 1:numSieve
%         results = [results; (sum(totpassInfo(k,:)))./sum(totVolume).*100];
%     end
%     semilogx(sieveSz,results,'-*')
%     ylim([0 100])
%     title('Grain Size Distribution')
%     xlabel('Grain Size Diameter (mm)')
%     ylabel('Cummulative (%) defined by Volume')
%     grid on
%     hold on
    catch
        disp('not a valid folder')
    end
end

figure(1)
results = (sum(totpassInfo,2)./sum(totVolume)) * 100;
semilogx(sieveSz,results,'-*')
ylim([0 100])
title('Grain Size Distribution')
xlabel('Grain Size (pixels)')
ylabel('Cummulative Passing (%), Volume')
grid on

function [pass, totVol] = processImg(filename,sieveSz)
BW=imread(filename); 
BW=imfill(BW,'holes');
BW=wiener2(BW,[5 5]);% you can change the size for  this one, e.g., [9,9]

% get properties of individual particles
props = regionprops('table',BW,'Centroid','MajorAxisLength',...
'MinorAxisLength','ConvexArea','Area','Perimeter');

%scale = .10769;
scale = 1;
major=table2array(props(:,{'MajorAxisLength'}))*scale;
minor=table2array(props(:,{'MinorAxisLength'}))*scale;
centroid=table2array(props(:,{'Centroid'}));
area=table2array(props(:,{'Area'}))*scale^2;
convex_area=table2array(props(:,{'ConvexArea'}));
perimeter=table2array(props(:,{'Perimeter'}));

grainsize = sqrt(0.5*(minor.^2 + major.^2));
volumeest = (4/3).*area.*(1/2).*(0.6837*minor);

%grain size vs. percent passing
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