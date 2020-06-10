function growAgg(aggRepo, cubeCell, scaleStepFactor)
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    scaleUp = [scaleStepFactor 0 0; 0 scaleStepFactor 0; 0 0 scaleStepFactor];
    scaleDown = [1/scaleStepFactor 0 0; 0 1/scaleStepFactor 0; 0 0 1/scaleStepFactor];
    
    for i = 1:numAgg
        curAggName = aggNames{i}
        curCubeNum = aggRepo.(curAggName).cubeNum;
        curCubePoints = cubeCell{curCubeNum};
        curCubeAlpha = alphaShape(curCubePoints(:,:,:));
        numTries = 1
        while 1
            aggRepo.(curAggName).Points = ...
                aggRepo.(curAggName).Points  * scaleUp
            pointCheck = ~inShape(curCubeAlpha, aggRepo.(curAggName).Points)
            if pointCheck > 0
                aggRepo.(curAggName).Points = ...
                    aggRepo.(curAggName).Points * scaleDown
                disp(numTries)
                break
            end
            numTries = numTries + 1;
        end
    end
end

