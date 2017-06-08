function nClustersOld = doPlots_TimeCourse(clusterIDs, data,featureNames, nClustersOld, sampleNumber, varargin)

uniqueClusters = unique(clusterIDs);
nClusters = length(uniqueClusters);


times = featureNames;
myCols = cool(nClusters);
nPlots = ceil(sqrt(nClusters));
if(size(varargin,2) == 0)
    nClustersOld = nClusters;
    clf;
    set(gcf,'NextPlot','add');
    axes;
%     h = title(['Iteration number: ', num2str(sampleNumber)]);
%     set(gca,'Visible','off');
%     set(h,'Visible','on');
else
    delete(findobj(varargin{1}, 'type', 'axes'))
    a = axes('parent',varargin{1});
%     h = title(a, ['Iteration number: ', num2str(sampleNumber)]);
%     set(a,'Visible','off');
%     set(h,'Visible','on');
% 
end
counter = 1;
hvec = zeros(1,length(uniqueClusters));
for i = uniqueClusters
    hvec(counter) = subplot(nPlots,nPlots,counter);
    plot(times, data(clusterIDs == i,:), 'color', myCols(counter,:));
    set(gca, 'Color', [.2 .2 .2],'YTickLabel',[],'XTickLabel',[],...
        'YTick',[],'XTick',[], 'xminortick','off', 'yminortick',...
        'off', 'xminorgrid','off', 'yminorgrid','off')
    xlim([min(times) max(times)]);
    counter = counter + 1;
end
set(get(hvec(1), 'Parent'), 'Title', ['Iteration number: ', num2str(sampleNumber)]);


end
