function [aggRepo, removeFields] = removeOverlaps(aggRepo)
%Removes overlaps between aggregates using alphaShape
%Input: A struct containing aggregates with overlaps
%Ouput: A struct containing aggregates with no overlaps

    %Setting constants
    aggN = fieldnames(aggRepo);
    numAggs = length(aggN);
    
    %Allocating cell to contain overlapping the smaller
    %overlapping aggregate
    removeFields = {};
    numRemove = 1;
    
    %iterating through each aggregate to see if there are any overlaps
    for i = 1:numAggs
            if any(strcmp(removeFields, aggN{i}))
                continue
            end
            curAggPoints = aggRepo.(aggN{i}).Points;
            
            %Getting the 20 closest aggregates
            aggNames = getClosestAggregates(getCentroid(curAggPoints), aggRepo);
            numCompare = length(aggNames);
           
         %iterating through the the closest aggregates and seeing if there
         %are any overlaps
        for j = 2:numCompare
            %Skipping any aggregates already removed
            if any(strcmp(removeFields, aggNames{j}))
                continue
            end
            otherAggPoints = aggRepo.(aggNames{j}).Points;
            crit = criticalAlpha(alphaShape(otherAggPoints), "one-region") + 30;
            otherAlpha = alphaShape(otherAggPoints, crit);
               
            %Using alphaShape's inShape function to detect overlaps
            %If there are any, the smaller aggregate is removed and added
            %to removeFields
            pointCheck = inShape(otherAlpha, curAggPoints);
            pointCheckSum = sum(pointCheck, 'all');
            if pointCheckSum > 0 %breaks if there is overlap
                curDiam = maxDiam(curAggPoints);
                otherDiam = maxDiam(otherAggPoints);
                if curDiam < otherDiam
                    removeFields{numRemove,1} = aggN{i};
                    numRemove = numRemove + 1;
                else
                    removeFields{numRemove,1} = aggNames{j};
                    numRemove = numRemove + 1;
                end
            end
        end
    end
    
    %remove aggregates from aggRepo
    aggRepo = rmfield(aggRepo, removeFields);
end

function closestAggs = getClosestAggregates(xyz, aggRepo)
    aggs = fieldnames(aggRepo);
    for i = 1:length(aggs)
        aggs(i,2) = mat2cell(norm(getCentroid(aggRepo.(aggs{i}).Points) - xyz), [1]);
    end
    aggs = sortrows(aggs, 2);
    closestAggs = aggs(1:20);
end 

function centroid = getCentroid(datapoints)
    x = datapoints(:,1);
    y = datapoints(:,2);
    z = datapoints(:,3);

    xcm = sum(x)./length(x);
    ycm = sum(y)./length(y);
    zcm = sum(z)./length(z);
    centroid = [xcm ycm zcm];
end