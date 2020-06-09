clear
clc

cubeSizeInit = 300;
nDivisions = 10;
step = cubeSizeInit/nDivisions;
nDivisions = nDivisions + 1;

cube = cell(nDivisions, nDivisions, nDivisions);
for x= 1:nDivisions
    for y = 1:nDivisions
        for z = 1:nDivisions
            cube{x,y,z} = [(step*(x-1)) (step*(y-1)) (step*(z-1))];
        end
    end
end

num = 1;
nDivisionsA = nDivisions-1
preCell = cell(nDivisionsA^3,1);
for x = 1:nDivisionsA
    for y = 1:nDivisionsA
        for z = 1:nDivisionsA
            coordCell = [cube{x,y,z}; cube{x,y+1,z};...
                                  cube{x+1,y,z}; cube{x+1,y+1,z};...
                                  cube{x,y,z+1}; cube{x,y+1,z+1};...
                                  cube{x+1,y,z+1}; cube{x+1,y+1,z+1}];
             preCell{num,1} = coordCell;
             num = num + 1;
        end
    end
end

cubeAlpha1 = alphaShape(preCell{1,1});
plot(cubeAlpha1);
hold on 
for i = 2:nDivisionsA^3
    plot(alphaShape(preCell{i,1}));
end
hold off