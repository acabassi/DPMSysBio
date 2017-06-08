function [finalClustering, PSM, reorderedPSM] = findSummaryClustering(fileName, uniqueIdentifier, burnin)
saveFileName = [strtok(fileName, '.'),'_Results_Chain', num2str(uniqueIdentifier)];

load(saveFileName)

allData = csvread([saveFileName '.csv'], 1,0);

alphas   = allData((burnin+1):end,1);
nSamples = length(alphas);
clusters = allData((burnin+1):end,2:(nGenes+1));

if(strcmp(dataType, 'TimeCourseUnequalSpacing'))
    hypers = allData((burnin+1):end,(nGenes+1):end);
end

% nGenes    = 400;
% nSamples  = 500;
% load('clsdraw')
% clusters  = clsdraw;
% 
allCombs  = combntns(1:nGenes,2);

clusters1 = clusters(:,allCombs(:,1));
clusters2 = clusters(:,allCombs(:,2));

PSM       = eye(nGenes);
PSM(sub2ind([nGenes nGenes], allCombs(:,1), allCombs(:,2))) = ...
    sum(clusters1 == clusters2)/nSamples;
sumpij    = sum(sum(PSM)) - nGenes;
PSM       = PSM + triu(PSM,1)';


histed    = histc(clusters', 1:max(max(clusters)));
myMat     = zeros(size(histed));


temp      = (arrayfun(@(x) nchoosek(x,2), histed(histed>=2))');
myMat     = zeros(size(histed));
myMat(histed>=2) = temp;
sumIij    = sum(myMat);


sumIpi = cellfun(@(x) compIpi(x,nGenes, PSM), num2cell(clusters,2));

correc = (sumIij.*sumpij)/nchoosek(nGenes,2);
pearsdraws = (sumIpi' - correc)./(0.5 * (sumpij + sumIij) - correc);


valdraws        = max(pearsdraws);
finalClustering = clusters(find(pearsdraws == valdraws,1),:);
[sorted, IX]    = sort(finalClustering);

reorderedPSM    = PSM(IX,IX);
figure
imagesc(reorderedPSM); colorbar
save([saveFileName 'summary.mat'], 'finalClustering', 'PSM', 'reorderedPSM');

doPlots(finalClustering, data,featureNames, nClustersOld, timeCourseSwitch, multinomialSwitch, bagOfWordsSwitch);
% uniqueClusterLabels = unique(finalClustering);
% for i = 1:length(uniqueClusterLabels)
%     figure
%     plot(
% end


end