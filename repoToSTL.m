function repoToSTL(aggRepo)
    folder = 'STL Files\Aggregates Out\'; %folder MUST be created
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    
    allPnts = aggRepo.(aggNames{1}).Points;
    allTris = aggRepo.(aggNames{1}).Faces;
    triLen = length(allPnts);
    for i = 2:numAgg %write to STL
        curAggName = aggNames{i};
        curAggPoints = aggRepo.(curAggName).Points;
        curAggFaces = aggRepo.(curAggName).Faces;
        
        addAggFaces = curAggFaces + triLen;
        allPnts = vertcat(allPnts, curAggPoints);
        allTris = vertcat(allTris, addAggFaces);
        
        triLen = triLen + length(curAggPoints);
        
        TR = triangulation(curAggFaces, curAggPoints);
        fn = strcat(folder, curAggName, ".stl");
        stlwrite(TR, fn);
    end
    
    allTR = triangulation(allTris, allPnts);
    allFn = strcat(folder, 'All', '.stl');
    stlwrite(allTR, allFn);
end