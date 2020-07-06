function plotSTL(folder)
    folder = strcat(folder, '\*.stl');
    fileNames = dir(folder);
    for i = 1:(length(fileNames))
        fn = fileNames(i,1).name;
        fn = strcat(folder(1:end-5), fn);
        aggSTL = stlread(fn);
        trimesh(aggSTL.ConnectivityList, aggSTL.Points(:,1), aggSTL.Points(:,2), aggSTL.Points(:,3))
        hold on
    end
end