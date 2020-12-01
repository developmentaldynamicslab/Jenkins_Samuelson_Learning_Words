function [ newCluster ] = rClusterTweak( leftCluster,rightCluster,comboHeight)
sizeL = size(leftCluster);
sizeR = size(rightCluster);
sizeL = sizeL(2);
sizeR = sizeR(2);
coin = randi(2);
if coin == 1
    rightCluster = -1*rightCluster;
else
    leftCluster = -1*leftCluster;
end
newCluster = [leftCluster rightCluster];
sizeNew = size(newCluster);
sizeNew = sizeNew(2);
avgLink = 0;
for l = 1:sizeNew
    if l ~= sizeNew
    for ll = l+1:sizeNew
        avgLink = avgLink + abs(newCluster(ll)-newCluster(l));
    end
    end
end
linksTotal = (sizeNew*(sizeNew-1))/2;
linksBetween = sizeL*sizeR;
avgLink = avgLink / linksTotal;
addend = (comboHeight-avgLink)/(linksBetween / linksTotal);
if coin == 1
    addend = -1*addend;
end
newRightCluster = rightCluster+addend;
minAll = min([leftCluster newRightCluster]);
newCluster = [leftCluster newRightCluster] - minAll + 1;

%{
%check
newCluster = [leftCluster newRightCluster];
sizeNew = size(newCluster);
sizeNew = sizeNew(2);
avgLink = 0;
for c = 1:sizeNew
    if c ~= sizeNew
    for cc = c+1:sizeNew
        avgLink = avgLink + abs(newCluster(cc)-newCluster(c));
    end
    end
end
avgLink = avgLink / linksTotal
%}
end

