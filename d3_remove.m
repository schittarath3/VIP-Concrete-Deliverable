clear
clc

%create psuedo binary image
mat1 = [1 1 0 0; 1 1 0 0; 0 0 0 0; 0 0 0 0];
mat2 = [0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 1 1];
mat3 = [0 0 0 0; 0 0 0 0; 0 0 1 1; 0 0 1 1];

%create 3D matrix
mat3D = cat(3, mat1, mat2, mat3);

CC = bwconncomp(mat3D, 6);