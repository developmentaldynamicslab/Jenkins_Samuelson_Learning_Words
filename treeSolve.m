function [solution] = treeSolve(adultdata, tree, n_field,constraint)
%treeSolve ...Put in a number for which tree matrix you want as input, and
%get a random, valid output in the single dimensional array format required
%by Bayes.  1 = animal 2 = vehicle 3 = vegetable
%constraint = 0 for all possibilities, 1 for subord on the left, 1 for
%subord in the middle.
reasonable = 0;
count = 0;
while (reasonable == 0) && count <10000
count = count + 1;
if adultdata == 1
if tree == 1
    %animals
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak([1,1],[1],0.011),[1],0.153),[1],0.170),rClusterTweak(1,1,0.153),0.193),1,0.324),1,0.336),1,0.364),1,0.418),rClusterTweak(1,1,0.397),0.443);
    solution = [raw(1) raw(1) raw(1) raw(4) raw(5) raw(11) raw(12) raw(2) raw(3) raw(6) raw(7) raw(8) raw(9) raw(10) raw(13)];
end
if tree == 2
    %vehicles
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.170),rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.04),1,0.067),1,0.153),0.191),1,0.247),1,0.278),1,0.353),rClusterTweak(1,rClusterTweak(rClusterTweak(1,1,0.170),1,0.250),0.378),0.4);
    solution = [raw(3) raw(3) raw(3) raw(1) raw(7) raw(10) raw(11) raw(4) raw(5) raw(2) raw(6) raw(8) raw(9) raw(12) raw(13)];
end
if tree == 3
    %vegetables
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.096),1,0.113),rClusterTweak(1,[1 1 1],0.042),0.116),1,0.289),rClusterTweak(rClusterTweak(1,1,0.290),rClusterTweak(1,1,0.246),0.314),0.344),1,0.368);
    solution = [raw(5) raw(5) raw(5) raw(1) raw(4) raw(9) raw(11) raw(6) raw(7) raw(2) raw(3) raw(8) raw(13) raw(12) raw(10)];
end
else

if tree == 1
    %animals
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak([1],[1],0.28),[1],0.281),[1],0.4325),[1],0.44),[1],0.475),[1],0.485),1,0.502),1,0.523),1,0.53),rClusterTweak(1,1,0.515),0.532),1,0.54);
    solution = [raw(1) raw(1) raw(1) raw(7) raw(8) raw(4) raw(5) raw(2) raw(3) raw(9) raw(10) raw(6) raw(11) raw(12) raw(13)];
end
if tree == 2
    %vehicles
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.215),1,0.235),rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.495),rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.297),1,0.457),1,0.5),0.501),rClusterTweak(rClusterTweak(1,1,0.385),1,0.449),0.51),0.523),1,0.585);
    solution = [raw(1) raw(1) raw(1) raw(4) raw(5) raw(6) raw(7) raw(2) raw(3) raw(10) raw(13) raw(8) raw(9) raw(11) raw(12)];
end
if tree == 3
    %vegetables
    raw = rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.2),rClusterTweak(1,1,0.22),0.442),rClusterTweak(rClusterTweak(rClusterTweak(1,1,0.2),1,0.256),1,0.395),0.48),1,0.515),1,0.53),rClusterTweak(1,1,0.323),0.582),1,0.62);
    solution = [raw(5) raw(5) raw(5) raw(1) raw(2) raw(8) raw(9) raw(6) raw(7) raw(3) raw(4) raw(10) raw(11) raw(12) raw(13)];
end
end
range = max(solution) - min(solution);
scalar = (n_field-40)/range;
solution = solution * scalar;
solution = solution - min(solution) + 20;


if constraint == 0 %don't check for reaosnable solutions.  otherwise, will run through and verify no weird overlaps, consistent with a groomed data set like in X&T
    reasonable = 1;
elseif constraint == 1
    if solution(8) < min(solution(10:15)) & solution(9) < min(solution(10:15)) & solution(10) < min(solution(12:15)) & solution(11) < min(solution(12:15)) & solution(1) < min(solution(4:5)) & solution(2) < min(solution(4:5)) & solution(4) < min(solution(6:7)) & solution(5) < min(solution(6:7))
        reasonable = 1;
    end
elseif constraint == 2
    if (abs(solution(10)-solution(8)) < abs(solution(12)-solution(8))) & (abs(solution(11)-solution(8)) < abs(solution(12)-solution(8))) & (abs(solution(10)-solution(8)) < abs(solution(13)-solution(8))) & (abs(solution(11)-solution(8)) < abs(solution(13)-solution(8))) & (abs(solution(10)-solution(8)) < abs(solution(14)-solution(8))) & (abs(solution(11)-solution(8)) < abs(solution(14)-solution(8))) & (abs(solution(10)-solution(8)) < abs(solution(15)-solution(8))) & (abs(solution(11)-solution(8)) < abs(solution(15)-solution(8))) & (abs(solution(4)-solution(1)) < abs(solution(6)-solution(1))) & (abs(solution(5)-solution(1)) < abs(solution(6)-solution(1))) & (abs(solution(4)-solution(1)) < abs(solution(7)-solution(1))) & (abs(solution(5)-solution(1)) < abs(solution(7)-solution(1)))
        reasonable = 1;
    end
else
    %figure out the largest distance between any two sub level items
    biggestSubDiff = 0;
    if (abs(solution(1)-solution(2)) > biggestSubDiff)
        biggestSubDiff = abs(solution(1)-solution(2));
    end
    if (abs(solution(1)-solution(3)) > biggestSubDiff)
        biggestSubDiff = abs(solution(1)-solution(3));
    end
    if (abs(solution(1)-solution(8)) > biggestSubDiff)
        biggestSubDiff = abs(solution(1)-solution(8));
    end
    if (abs(solution(1)-solution(9)) > biggestSubDiff)
        biggestSubDiff = abs(solution(1)-solution(9));
    end
    if (abs(solution(2)-solution(3)) > biggestSubDiff)
        biggestSubDiff = abs(solution(2)-solution(3));
    end
    if (abs(solution(2)-solution(8)) > biggestSubDiff)
        biggestSubDiff = abs(solution(2)-solution(8));
    end
    if (abs(solution(2)-solution(9)) > biggestSubDiff)
        biggestSubDiff = abs(solution(2)-solution(9));
    end
    if (abs(solution(3)-solution(8)) > biggestSubDiff)
        biggestSubDiff = abs(solution(3)-solution(8));
    end
    if (abs(solution(3)-solution(9)) > biggestSubDiff)
        biggestSubDiff = abs(solution(3)-solution(9));
    end
    if (abs(solution(8)-solution(9)) > biggestSubDiff)
        biggestSubDiff = abs(solution(8)-solution(9));
    end
    
    %check each basic level or higher item to see if any are closer than
    %this to any sub items
    tooNears = 0;
    if (abs(solution(1) - solution(4)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(5)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(10)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(11)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(6)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(7)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(12)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(13)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(14)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(1) - solution(15)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(4)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(5)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(10)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(11)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(6)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(7)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(12)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(13)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(14)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(2) - solution(15)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(4)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(5)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(10)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(11)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(6)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(7)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(12)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(13)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(14)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(3) - solution(15)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(4)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(5)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(10)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(11)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(6)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(7)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(12)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(13)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(14)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(8) - solution(15)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(4)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(5)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(10)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(11)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(6)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(7)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(12)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(13)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(14)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    if (abs(solution(9) - solution(15)) < biggestSubDiff)
        tooNears = tooNears + 1;
    end
    %if so, count how MANY sub items they are closer to than this at once.
    %if > 2, not reasonable .
    %tooNears
    %count
    if (tree == 1 && tooNears < 8) %8 was used during Omega tweaking, 3 minimum, will always give same results if chosen
        tooNears;
        reasonable = 1;
    end
    if (tree == 2 && tooNears < 17) %17 was used during Omega tweaking, 12 minimum
        tooNears;
        reasonable = 1;
    end
    if (tree == 3 && tooNears < 10) %10 was used during Omega tweaking, 5 minimum
        tooNears;
        reasonable = 1;
    end
end

%{
swapM=zeros(1,5);
swapM = randswap([solution(1) solution(2) solution(3) solution(8) solution(9)]);
solution(1)=swapM(1);
solution(2)=swapM(2);
solution(3)=swapM(3);
solution(8)=swapM(4);
solution(9)=swapM(5);

swapM=zeros(1,4);
swapM = randswap([solution(4) solution(5) solution(10) solution(11)]);
solution(4)=swapM(1);
solution(5)=swapM(2);
solution(10)=swapM(3);
solution(11)=swapM(4);

swapM=zeros(1,6);
swapM = randswap([solution(6) solution(7) solution(12) solution(13) solution(14) solution(15)]);
solution(6)=swapM(1);
solution(7)=swapM(2);
solution(12)=swapM(3);
solution(13)=swapM(4);
solution(14)=swapM(5);
solution(15)=swapM(6);
%}

end %end of the while loop



