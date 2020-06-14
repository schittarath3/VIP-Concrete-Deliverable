function volFract = volFracCheck(aggRepo, cubeLength)
%Gets volume fraction of aggregate within the cube
%Inputs: 
%   aggRepo: Struct of aggregates of form aggName ->cubeNum,Faces,Points
%   cubeLength: Length of cube
%Outputs:
%   Volume fraction of aggregates within cube as a float
    cubeVol = cubeLength^3; %Volume of cube
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    totalAggVol = 0;
    for i = 1:numAgg %Adds the volume of aggregate alphashapes
        curAgg = aggRepo.(aggNames{i}).Points;
        aggAlpha = alphaShape(curAgg);
        aggVol = volume(aggAlpha);
        totalAggVol = totalAggVol + aggVol;
    end
    volFract = totalAggVol/cubeVol;
end