function repoToSTL(aggRepo, iter)
    folder = 'STL Files\Aggregates Out\'; %folder MUST be created
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
     if iter == 1
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
     else
        allPnts = aggRepo.(aggNames{1}).OriginalPoints;
        allTris = aggRepo.(aggNames{1}).OriginalFaces;
        triLen = length(allPnts);
        for i = 2:numAgg %write to STL
            curAggName = aggNames{i};
            curAggPoints = aggRepo.(curAggName).OriginalPoints;
            curAggFaces = aggRepo.(curAggName).OriginalFaces;

            addAggFaces = curAggFaces + triLen;
            allPnts = vertcat(allPnts, curAggPoints);
            allTris = vertcat(allTris, addAggFaces);

            triLen = triLen + length(curAggPoints);

            TR = triangulation(curAggFaces, curAggPoints);
            fn = strcat(folder, curAggName, ".stl");
            stlwrite(TR, fn);
        end
     end

    allTR = triangulation(allTris, allPnts);
    allFn = strcat(folder, 'All', '.stl');
    stlwrite(allTR, allFn);
end