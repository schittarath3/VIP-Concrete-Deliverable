function dist = maxDiam(aggpts)
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