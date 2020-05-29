%comment
fileList = dir('STL Files\Aggregates Processed\*.stl');
fileList = struct2cell(fileList);
aggList = cell(1, length(fileList));
samples = length(fileList);
meshvol = zeros(1,samples);
for i = 1:samples
    %obtain stl
    folder = fileList{2, i};
    file = fileList{1, i};
    fn = strcat(folder, '\', file);
    aggList{i} = stlread(fn);
end