function nClustersOld = doPlots_BagOfWords(clusterIDs, data,featureNames, nClustersOld, sampleNumber, varargin)
if(size(varargin,2) ~= 0)
    delete(findobj(varargin{1}, 'type', 'axes'))
    a = axes('parent',varargin{1});
end

uniqueClusters = unique(clusterIDs);
nClusters = length(uniqueClusters);




mygray = gray;
mygray = mygray(end:-1:1,:);




counter = 1;
newData = [];
for i = uniqueClusters
    tempData = data(clusterIDs == i,:);
    if(mod(counter,2))
        tempData(tempData == 0) = 0.15;
    end
    newData  = [newData; tempData];
    counter = counter + 1;
end

imagesc(newData')
colormap(mygray);

tickPoints = cumsum(histc(clusterIDs, uniqueClusters));

set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[])
set(gca,'XTick',tickPoints+0.5, 'TickLength', [.4 .4],'LineWidth', 1)
if(size(varargin,2) == 0)
    title(['Iteration number: ', num2str(sampleNumber)])
else
    set(varargin{1}, 'Title', ['Iteration number: ', num2str(sampleNumber)]);
end


end

