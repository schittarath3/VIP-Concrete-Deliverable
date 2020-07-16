% Grain size distribution from 2D binary image of rocks using watershed
% segmentation algorithm
% You will need image processing license to run this code
% 
% Please cite this paper: 
% Rabbani, A., Ayatollahi, S. (2015). Comparing three image processing algorithms 
% to estimate the grain-size distribution of porous rocks from binary 2d images and 
% sensitivity analysis of the grain overlapping degree. Special Topics & Reviews in 
% Porous Media: An International Journal 6 (1), 71-89.
% 
% By Arash Rabbani
% arashrabbani.com
% rabarash@yahoo.com
clc;clear; close all;
% INPUTS
A=imread('bl_top0054_post.png');
Resolution=5; % micron/pixel % this is the spatial resolution of the input 
Bins=20; % This is the number of bars for histogram chart
% CALCULATIONS
Conn=8;
[s1,s2]=size(A);
A=~bwmorph(A,'majority',10);
Poro=sum(sum(~A))/(s1*s2);
D=-bwdist(A,'cityblock');
B=medfilt2(D,[3 3]);
B=watershed(B,Conn);
Pr=zeros(s1,s2);
for I=1:s1
    for J=1:s2
        if A(I,J)==0 && B(I,J)~=0
            Pr(I,J)=1;
        end
    end
end
Pr=bwareaopen(Pr,9,Conn);
[Pr_L,Pr_n]=bwlabel(Pr,Conn);
V=zeros(Pr_n,1);
for I=1:s1
    for J=1:s2
        if Pr_L(I,J)~=0
            V(Pr_L(I,J))=V(Pr_L(I,J))+1;
        end
    end
end
R=Resolution.*(V./pi).^.5; % grain radius
%Outputs
Average_grain_radius_micron=mean(R)
Standard_deviation_of_grain_radius_micron=std(R)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
RGB=label2rgb(Pr_L,'jet', 'w', 'shuffle');
imshow(RGB)
imwrite(RGB,'Output.png')
subplot(1,2,2)
Rel_Frequencies=hist(R,[1:round(max(R)/Bins):round(max(R))])./sum(sum(hist(R,[1:round(max(R)/Bins):round(max(R))]))); 
bar([1:round(max(R)/Bins):round(max(R))],Rel_Frequencies); 
xlabel('Equivalent Grain Radius (micron)'); ylabel('Relative Frequency'); axis([1 max(R) 0 max(Rel_Frequencies)]); axis square;
annotation('textbox',[.2 .85 .1 .1], 'String', [ 'Average grain radius = ' num2str(Average_grain_radius_micron) ' micron'])