function nClustersOld = doPlots_Multinomial(clusterIDs, data,featureNames, nClustersOld, sampleNumber, varargin)


if(size(varargin,2) ~= 0)
    delete(findobj(varargin{1}, 'type', 'axes'))
    a = axes('parent',varargin{1});
end
uniqueClusters = unique(clusterIDs);


newData = [];
for i = uniqueClusters
    newData = [newData; data(clusterIDs == i,:)];
end
h = imagesc(newData');
if(length(unique(newData(:))) == 3)
    mycolormap = [0 0 1; 1 1 1; 1 0 0];
elseif(length(unique(newData(:))) == 2)
    mycolormap  = [1 1 1; 0 0 0];
else
    mycolormap = jet(length(unique(newData(:))));
end
colormap(mycolormap);
tickPoints = cumsum(histc(clusterIDs, uniqueClusters));

set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[])
set(gca,'XTick',tickPoints+0.5, 'TickLength', [.5 .5],'LineWidth', 1.5)

% if(length(unique(newData(:))) == 3)
%     colorbar('YTick', [1.333 2 2.66], 'YTickLabel', [1, 2, 3], 'TickLength',[ 0 0 ])
% else
maxNum = max(newData(:)); minNum = min(newData(:));
minmaxDiff = maxNum - minNum;
stepsize   = minmaxDiff/maxNum;
yticklabels = minNum:maxNum;
yticks = (minNum+stepsize/2):stepsize:maxNum;
colorbar('YTick', yticks, 'YTickLabel', yticklabels, 'TickLength', 0)
% end

if(size(varargin,2) == 0)
    title(['Iteration number: ', num2str(sampleNumber)])
else
    set(varargin{1}, 'Title', ['Iteration number: ', num2str(sampleNumber)]);
end

end
