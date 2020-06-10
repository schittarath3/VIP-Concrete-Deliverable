function newRepo = growAgg(aggRepo, cubeCell, scaleStepFactor)
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    scaleUp = [scaleStepFactor 0 0; 0 scaleStepFactor 0; 0 0 scaleStepFactor];
    scaleDown = [1/scaleStepFactor 0 0; 0 1/scaleStepFactor 0; 0 0 1/scaleStepFactor];
    
    newRepo = struct
    for i = 1:numAgg
        curAggName = aggNames{i};
        curCubeNum = aggRepo.(curAggName).cubeNum;
        curCubePoints = cubeCell{curCubeNum};
        curCubeAlpha = alphaShape(curCubePoints(:,:,:));
        numTries = 1;
        newRepo.(curAggName).cubeNum = curCubeNum;
        newRepo.(curAggName).Points = aggRepo.(curAggName).Points
        while 1
            newRepo.(curAggName).Points = ...
                newRepo.(curAggName).Points  * scaleUp;
            pointCheck = ~inShape(curCubeAlpha, newRepo.(curAggName).Points);
            if pointCheck > 0
                newRepo.(curAggName).Points = ...
                    newRepo.(curAggName).Points * scaleDown;
                disp(numTries)
                break
            end
            numTries = numTries + 1;
        end
    end
end

