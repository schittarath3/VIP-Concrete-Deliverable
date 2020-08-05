function dist = maxDiam(aggpts)
%Find the longest length of an aggregate by selecting a fixed point and
%finding the distance using the vector norm. The distance is compared to
%the original distance in each iteration, updating the distance vector with
%the longest length.
%Inputs: 
%   aggpts - matrix of aggregates coordinates with each column representing
%   (x,y,z) coordinates (generated from reading stl file of the aggregate).
%Outputs:
%   dist - the maximum length of the aggregates.

dist = 0;
for i = 1:length(aggpts)
    for k = 1:length(aggpts)
    distp = abs(norm(aggpts(i,:) - aggpts(k,:)));
    if distp > dist
        dist = distp;
    end
    end
end
end