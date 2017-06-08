
function nClustersOld = doPlots(clusterIDs, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch)

uniqueClusters = unique(clusterIDs);
nClusters = length(uniqueClusters);


if(timeCourseSwitch)
    times = featureNames;
    myCols = cool(nClusters);
    nPlots = ceil(sqrt(nClusters));
    if(nClusters ~= nClustersOld)
        clf;
        nClustersOld = nClusters;
    end
    counter = 1;
    for i = uniqueClusters
        subplot(nPlots,nPlots,counter)
        plot(times, data(clusterIDs == i,:), 'color', myCols(counter,:));
        set(gca, 'Color', [.3 .3 .3],'YTickLabel',[],'XTickLabel',[],...
            'YTick',[],'XTick',[], 'xminortick','off', 'yminortick',...
            'off', 'xminorgrid','off', 'yminorgrid','off')
        xlim([min(times) max(times)]);
        counter = counter + 1;
    end
end
if(multinomialSwitch)
    newData = [];
    for i = uniqueClusters
        newData = [newData; data(clusterIDs == i,:)];
    end
    imagesc(newData');
    if(length(unique(newData(:))) == 3)
        mycolormap = [0 0 1; 1 1 1; 1 0 0];
        colormap(mycolormap);
    else
        colormap jet;
    end
    tickPoints = cumsum(histc(clusterIDs, uniqueClusters));
    
    set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[])
    set(gca,'XTick',tickPoints+0.5, 'TickLength', [.4 .4],'LineWidth', 1.5)
end
if(bagOfWordsSwitch)
    
    mygray = gray;
    mygray = mygray(end:-1:1,:);
    
    
    
    
    counter = 1;
    allClusterData = [];
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
    
    
    
end
end
