function repoToSTL(aggRepo)
    folder = 'STL Files\Aggregates Out\'; %folder MUST be created
    aggNames = fieldnames(aggRepo);
    numAgg = length(aggNames);
    for i = 1:numAgg %write to STL
        curAggName = aggNames{i};
        curAggPoints = aggRepo.(curAggName).Points;
        curAggFaces = aggRepo.(curAggName).Faces;
        
        TR = triangulation(curAggFaces, curAggPoints);
        fn = strcat(folder, curAggName, ".stl");
        stlwrite(TR, fn);
    end
end