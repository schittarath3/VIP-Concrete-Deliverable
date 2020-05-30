%comment
fileList = dir('STL Files\Aggregates Processed 3\*.stl')
fileList = struct2cell(fileList);
aggList = cell(1, length(fileList));
samples = size(fileList);
meshvol = zeros(1,samples(2));
for i = 1:samples(2)
    %obtain stl
    folder = fileList{2, i};
    file = fileList{1, i};
    fn = strcat(folder, '\', file);
    aggList{i} = stlread(fn);
end